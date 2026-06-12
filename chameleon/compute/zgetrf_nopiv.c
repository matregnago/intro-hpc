/**
 *
 * @file zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv wrappers
 *
 * @version 1.4.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Pierre Esterie
 * @date 2025-12-19
 *
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Allocate the required workspaces for asynchronous getrf
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 *******************************************************************************
 *
 * @retval An allocated opaque pointer to use in CHAMELEON_zgetrf_nopiv_Tile_Async()
 * and to free with CHAMELEON_zgetrf_nopiv_WS_Free().
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_zgetrf_nopiv_WS_Free
 *
 */
void *CHAMELEON_zgetrf_nopiv_WS_Alloc( const CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    struct chameleon_pzgetrf_nopiv_s *options;
    int P, Q;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return NULL;
    }

    options = calloc( 1, sizeof(struct chameleon_pzgetrf_nopiv_s) );
    options->use_workspace = 0;

    P = chameleon_desc_datadist_get_iparam( A, 0 );
    Q = chameleon_desc_datadist_get_iparam( A, 1 );

    if ( ( ( P > 1 ) || ( Q > 1 ) ) &&
         ( A->get_rankof_init == chameleon_getrankof_2d ) &&
         ( chamctxt->generic_enabled != CHAMELEON_TRUE ) )
    {
        int lookahead = chamctxt->lookahead;
        options->use_workspace = 1;

        chameleon_desc_init( chamctxt, &(options->WL), "GETRF_NP_WL", CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, A->mb, A->nb,
                             A->mt * A->mb, A->nb * Q * lookahead,
                             A->mt * A->mb, A->nb * Q * lookahead, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );

        chameleon_desc_init( chamctxt, &(options->WU), "GETRF_NP_WU", CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, A->mb, A->nb,
                             A->mb * P * lookahead, A->nt * A->nb,
                             A->mb * P * lookahead, A->nt * A->nb, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );
    }

    return (void*)options;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Free the allocated workspaces for asynchronous getrf
 *
 *******************************************************************************
 *
 * @param[in,out] user_ws
 *          On entry, the opaque pointer allocated by CHAMELEON_zgetrf_nopiv_WS_Alloc()
 *          On exit, all data are freed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_zgetrf_nopiv_WS_Alloc
 *
 */
void CHAMELEON_zgetrf_nopiv_WS_Free( void *user_ws )
{
    struct chameleon_pzgetrf_nopiv_s *ws = (struct chameleon_pzgetrf_nopiv_s*)user_ws;

    if ( ws->use_workspace ) {
        chameleon_desc_destroy( &(ws->WL) );
        chameleon_desc_destroy( &(ws->WU) );
    }
    free( ws );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgetrf_nopiv - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm without row pivoting.
 *  WARNING: Don't use this function if you are not sure your matrix is diagonal
 *  dominant.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been
 *               completed, but the factor U is exactly singular, and division
 *               by zero will occur if it is used to solve a system of
 *               equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_nopiv_Tile
 * @sa CHAMELEON_zgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_cgetrf_nopiv
 * @sa CHAMELEON_dgetrf_nopiv
 * @sa CHAMELEON_sgetrf_nopiv
 *
 */
int CHAMELEON_zgetrf_nopiv( int M, int N,
                            CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    void *ws = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_nopiv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgetrf_nopiv", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgetrf_nopiv", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgetrf_nopiv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGESV, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrf_nopiv", "chameleon_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( &descAt );
    CHAMELEON_zgetrf_nopiv_Tile_Async( &descAt, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    CHAMELEON_zgetrf_nopiv_WS_Free( ws );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgetrf_nopiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of CHAMELEON_zgetrf_nopiv().  Operates on matrices stored by
 *  tiles.  All matrices are passed through descriptors.  All dimensions are
 *  taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_nopiv
 * @sa CHAMELEON_zgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_cgetrf_nopiv_Tile
 * @sa CHAMELEON_dgetrf_nopiv_Tile
 * @sa CHAMELEON_sgetrf_nopiv_Tile
 * @sa CHAMELEON_zgetrs_Tile
 *
 */
int CHAMELEON_zgetrf_nopiv_Tile( CHAM_desc_t *A )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    void               *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( A );
    CHAMELEON_zgetrf_nopiv_Tile_Async( A, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zgetrf_nopiv_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zgetrf_nopiv_Tile_Async - Computes the tile LU factorization of a
 *  matrix.  Non-blocking equivalent of CHAMELEON_zgetrf_nopiv_Tile().  May return
 *  before the computation is finished.  Allows for pipelining of operations ar
 *  runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_nopiv
 * @sa CHAMELEON_zgetrf_nopiv_Tile
 * @sa CHAMELEON_cgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_dgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_sgetrf_nopiv_Tile_Async
 * @sa CHAMELEON_zgetrs_Tile_Async
 *
 */
int CHAMELEON_zgetrf_nopiv_Tile_Async( CHAM_desc_t        *A,
                                       void               *user_ws,
                                       RUNTIME_sequence_t *sequence,
                                       RUNTIME_request_t  *request )
{
    CHAM_context_t                   *chamctxt;
    struct chameleon_pzgetrf_nopiv_s *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_nopiv_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_nopiv_Tile", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrf_nopiv_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zgetrf_nopiv_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

   if ( user_ws == NULL ) {
        ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( A );
    }
    else {
        ws = user_ws;
    }

    chameleon_pzgetrf_nopiv( ws, A, sequence, request );

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zgetrf_nopiv_WS_Free( ws );
    }
    return CHAMELEON_SUCCESS;
}
