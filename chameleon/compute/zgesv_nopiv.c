/**
 *
 * @file zgesv_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesv_nopiv wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgesv_nopiv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *  The tile LU decomposition with partial tile pivoting and row interchanges is used to factor A.
 *  The factored form of A is then used to solve the system of equations A * X = B.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] descL
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          necessary to solve the system of equations.
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesv_nopiv_Tile
 * @sa CHAMELEON_zgesv_nopiv_Tile_Async
 * @sa CHAMELEON_cgesv_nopiv
 * @sa CHAMELEON_dgesv_nopiv
 * @sa CHAMELEON_sgesv_nopiv
 *
 */
int CHAMELEON_zgesv_nopiv( int N, int NRHS,
                           CHAMELEON_Complex64_t *A, int LDA,
                           CHAMELEON_Complex64_t *B, int LDB )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    void *ws = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "illegal value of NRHS");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "illegal value of LDA");
        return -4;
    }
    if (LDB < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGESV, N, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgesv_nopiv", "chameleon_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( &descAt );
    CHAMELEON_zgesv_nopiv_Tile_Async( &descAt, &descBt, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    CHAMELEON_zgetrf_nopiv_WS_Free( ws );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgesv_nopiv_Tile - Solves a system of linear equations using the tile LU factorization.
 *  Tile equivalent of CHAMELEON_zgetrf_nopiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesv_nopiv
 * @sa CHAMELEON_zgesv_nopiv_Tile_Async
 * @sa CHAMELEON_cgesv_nopiv_Tile
 * @sa CHAMELEON_dgesv_nopiv_Tile
 * @sa CHAMELEON_sgesv_nopiv_Tile
 * @sa CHAMELEON_zcgesv_Tile
 *
 */
int CHAMELEON_zgesv_nopiv_Tile( CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    void               *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgesv_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( A );
    CHAMELEON_zgesv_nopiv_Tile_Async( A, B, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

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
 *  CHAMELEON_zgesv_nopiv_Tile_Async - Solves a system of linear equations using the tile
 *  LU factorization.
 *  Non-blocking equivalent of CHAMELEON_zgesv_nopiv_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
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
 * @sa CHAMELEON_zgesv_nopiv
 * @sa CHAMELEON_zgesv_nopiv_Tile
 * @sa CHAMELEON_cgesv_nopiv_Tile_Async
 * @sa CHAMELEON_dgesv_nopiv_Tile_Async
 * @sa CHAMELEON_sgesv_nopiv_Tile_Async
 * @sa CHAMELEON_zcgesv_Tile_Async
 *
 */
int CHAMELEON_zgesv_nopiv_Tile_Async( CHAM_desc_t        *A,
                                      CHAM_desc_t        *B,
                                      void               *user_ws,
                                      RUNTIME_sequence_t *sequence,
                                      RUNTIME_request_t *request )
{
    CHAM_context_t                   *chamctxt;
    struct chameleon_pzgetrf_nopiv_s *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgesv_nopiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgesv_nopiv_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgesv_nopiv_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgesv_nopiv_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgesv_nopiv_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        chameleon_error("CHAMELEON_zgesv_nopiv_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return CHAMELEON_SUCCESS;
     */

   if ( user_ws == NULL ) {
        ws = CHAMELEON_zgetrf_nopiv_WS_Alloc( A );
    }
    else {
        ws = user_ws;
    }
    chameleon_pzgetrf_nopiv( ws, A, sequence, request );

    chameleon_pztrsm( ChamLeft, ChamLower, ChamNoTrans, ChamUnit, (CHAMELEON_Complex64_t)1.0, A, B, sequence, request );

    chameleon_pztrsm( ChamLeft, ChamUpper, ChamNoTrans, ChamNonUnit, (CHAMELEON_Complex64_t)1.0, A, B, sequence, request );

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( B, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zgetrf_nopiv_WS_Free( ws );
    }
    return CHAMELEON_SUCCESS;
}
