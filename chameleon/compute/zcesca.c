/**
 *
 * @file zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca wrappers
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> s d c z
 *
 */
#include "control/common.h"
#include <stdlib.h>

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zcesca_WS_Alloc - Allocate the required workspaces for asynchronous center
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 *******************************************************************************
 *
 * @retval An allocated opaque pointer to use in CHAMELEON_zcesca_Tile_Async()
 * and to free with CHAMELEON_zcesca_WS_Free().
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zcesca_Tile_Async
 * @sa CHAMELEON_zcesca_WS_Free
 *
 */
void *CHAMELEON_zcesca_WS_Alloc( const CHAM_desc_t *A )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzcesca_s *options;
    int                         workmt, worknt;
    int                         P, Q;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return NULL;
    }

    options = calloc( 1, sizeof(struct chameleon_pzcesca_s) );

    P = chameleon_desc_datadist_get_iparam(A, 0);
    Q = chameleon_desc_datadist_get_iparam(A, 1);
    workmt = chameleon_max( A->mt, P );
    worknt = chameleon_max( A->nt, Q );

    chameleon_desc_init_2dtile( chamctxt, &(options->Wgcol), "CESCA_Wgcol", ChamComplexDouble,
                                1, A->nb, workmt, A->n, P, Q );

    chameleon_desc_init_2dtile( chamctxt, &(options->Wgrow), "CESCA_Wgrow", ChamComplexDouble,
                                A->mb, 1, A->m, worknt, P, Q );

    chameleon_desc_init_2dtile( chamctxt, &(options->Wgelt), "CESCA_Wgelt", ChamComplexDouble,
                                1, 1, 1, worknt, P, Q );

    chameleon_desc_init_2dtile( chamctxt, &(options->Wdcol), "CESCA_Wdcol", ChamRealDouble,
                                2, A->nb, 2*workmt, A->n, P, Q );

    chameleon_desc_init_2dtile( chamctxt, &(options->Wdrow), "CESCA_Wdrow", ChamRealDouble,
                                A->mb, 2, A->m, 2*worknt, P, Q );

    return (void*)options;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Free the allocated workspaces for asynchronous center
 *
 *******************************************************************************
 *
 * @param[in,out] user_ws
 *          On entry, the opaque pointer allocated by CHAMELEON_zcesca_WS_Alloc()
 *          On exit, all data are freed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zcesca_Tile_Async
 * @sa CHAMELEON_zcesca_WS_Alloc
 *
 */
void CHAMELEON_zcesca_WS_Free( void *user_ws )
{
    struct chameleon_pzcesca_s *ws = (struct chameleon_pzcesca_s*)user_ws;

    chameleon_desc_destroy( &(ws->Wgcol) );
    chameleon_desc_destroy( &(ws->Wgrow) );
    chameleon_desc_destroy( &(ws->Wgelt) );
    chameleon_desc_destroy( &(ws->Wdcol) );
    chameleon_desc_destroy( &(ws->Wdrow) );
    free( ws );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zcesca replace a general matrix by the Centered-Scaled matrix inplace.
 *  This algorithm is used as a pretreatment of a Principal Component Algorithm
 *  (PCA) or a Correspondence analysis (COA).
 *
 *  Considering a matrix A of size m x n, \f[A = (a_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f]
 *  Lets
 *  \f[g_i = \frac{1}{n} \sum_j a_{ij} \\
 *     g_j = \frac{1}{m} \sum_i a_{ij} \\
 *     g   = \frac{1}{mn} \sum_{i,j} a_{ij}\f]
 *  A centered rowwise gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_i \f]
 *  A centered columnwise gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_j \f]
 *  A bicentered gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_i - g_j + g \f]
 * Lets
 * \f[d_i = || a_{i*} || = \sqrt{ \sum_j a_{ij}^2} \\
 *    d_j = || a_{*j} || = \sqrt{ \sum_i a_{ij}^2} \f]
 * A scaled rowwise gives \f[A' = (a_{i,j}')_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ a_{i*}' = \frac{a_{i*}}{d_i} \f]
 * A scaled columnwise gives \f[A' = (a_{i,j}')_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ a_{*j}' = \frac{a_{*j}}{d_j} \f]
 *
 * This function can also be used to compute a pretreatment of a Correspondence analysis (COA).
 * To use it set center = 1, scale = 1, axis = ChamEltwise.
 * A on entry is a contingency table. For this pre-treatment we need to work on
 * the frequencies table (each value of A is divided by the global sum)
 * In this case lets
 * \f[r_i = \sum_j a_{ij} \\
 *    c_j = \sum_i a_{ij} \\
 *    sg  = \sum_{i,j} a_{ij} \f]
 *
 * A transformed gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ \bar{a}_{i,j} = \frac{a_{i,j}-r_i*c_j/sg}{ \sqrt{r_i*c_j} } \f]
 *
 * It also gives \f[ r_i \f] and \f[ c_j \f]
 * in vectors SR and SC useful in the post-treatment of the COA.
 * SR and SC must be already allocated.
 *
 *******************************************************************************
 *
 * @param[in] center
 *          1 if A must be centered, else 0.
 *
 * @param[in] scale
 *          1 if A must be scaled, else 0.
 *
 * @param[in] axis
 *          Specifies the axis over which to center and or scale.
 *            = ChamColumnwise: centered column-wise
 *            = ChamRowwise: centered row-wise
 *            = ChamEltwise:
 *              bi-centered if center=1 and scale=0
 *              COA if center=1 and scale=1
 *
 * @param[in] M
 *          The number of rows of the overall matrix.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the overall matrix.  N >= 0.
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] SR
 *          The vector of size M containing the \f[ r_i \f]
 *
 * @param[out] SC
 *          The vector of size N containing the \f[ c_j \f]
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zcesca_Tile
 * @sa CHAMELEON_zcesca_Tile_Async
 * @sa CHAMELEON_scesca
 *
 */
int CHAMELEON_zcesca(int center, int scale, cham_store_t axis,
                     int M, int N, CHAMELEON_Complex64_t *A, int LDA,
                     CHAMELEON_Complex64_t *SR, CHAMELEON_Complex64_t *SC)
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    struct chameleon_pzcesca_s *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zcesca", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (center != 0) && (center != 1) ) {
        chameleon_error("CHAMELEON_zcesca", "Illegal value of center");
        return -1;
    }
    if ( (scale != 0) && (scale != 1) ) {
        chameleon_error("CHAMELEON_zcesca", "Illegal value of scale");
        return -2;
    }
    if ( (axis != ChamColumnwise) && (axis != ChamRowwise) && (axis != ChamEltwise) ) {
        chameleon_error("CHAMELEON_zcesca", "Illegal value of axis");
        return -3;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zcesca", "Illegal value of M");
        return -4;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zcesca", "Illegal value of N");
        return -5;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zcesca", "illegal value of LDA");
        return -7;
    }

    /* Quick return */
    if (N == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_DGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zcesca", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    ws = CHAMELEON_zcesca_WS_Alloc( &descAt );
    CHAMELEON_zcesca_Tile_Async( center, scale, axis, &descAt, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* pre-coa case : save the sums over rows and columns */
    if ( (center == 1) && (scale == 1) && (axis == ChamEltwise) ) {
        CHAM_desc_t *descSR = chameleon_desc_submatrix( &(ws->Wgrow), 0, 0, M, 1 );
        CHAM_desc_t *descSC = chameleon_desc_submatrix( &(ws->Wgcol), 0, 0, 1, N );
        CHAMELEON_zDesc2Lap( ChamUpperLower, descSR, SR, M );
        CHAMELEON_zDesc2Lap( ChamUpperLower, descSC, SC, 1 );
        free( descSR );
        free( descSC );
    }

    /* Cleanup the temporary data */
    CHAMELEON_zcesca_WS_Free( ws );
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
 *  CHAMELEON_zcesca_Tile - Tile equivalent of CHAMELEON_zcesca().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] center
 *          1 if A must be centered, else 0.
 *
 * @param[in] scale
 *          1 if A must be scaled, else 0.
 *
 * @param[in] axis
 *          Specifies the axis over which to center and or scale.
 *            = ChamColumnwise: centered column-wise
 *            = ChamRowwise: centered row-wise
 *            = ChamEltwise:
 *              bi-centered if center=1 and scale=0
 *              pre-coa if center=1 and scale=1
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[out] SR
 *          The vector of size M containing the \f[ r_i \f]
 *
 * @param[out] SC
 *          The vector of size N containing the \f[ c_j \f]
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zcesca
 * @sa CHAMELEON_zcesca_Tile_Async
 * @sa CHAMELEON_scesca_Tile
 *
 */
int CHAMELEON_zcesca_Tile( int center, int scale, cham_store_t axis, CHAM_desc_t *A,
                           CHAMELEON_Complex64_t *SR, CHAMELEON_Complex64_t *SC )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;
    struct chameleon_pzcesca_s *ws;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zcesca_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (center != 0) && (center != 1) ) {
        chameleon_error("CHAMELEON_zcesca_Tile", "Illegal value of center");
        return -1;
    }
    if ( (scale != 0) && (scale != 1) ) {
        chameleon_error("CHAMELEON_zcesca_Tile", "Illegal value of scale");
        return -2;
    }
    if ( (axis != ChamColumnwise) && (axis != ChamRowwise) && (axis != ChamEltwise) ) {
        chameleon_error("CHAMELEON_zcesca_Tile", "Illegal value of axis");
        return -3;
    }

    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zcesca_WS_Alloc( A );
    CHAMELEON_zcesca_Tile_Async( center, scale, axis, A, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* pre-coa case : save the sums over rows and columns */
    if ( (center == 1) && (scale == 1) && (axis == ChamEltwise) ) {
        CHAM_desc_t *descSR = chameleon_desc_submatrix( &(ws->Wgrow), 0, 0, A->lm, 1 );
        CHAM_desc_t *descSC = chameleon_desc_submatrix( &(ws->Wgcol), 0, 0, 1, A->ln );
        CHAMELEON_zDesc2Lap( ChamUpperLower, descSR, SR, A->lm );
        CHAMELEON_zDesc2Lap( ChamUpperLower, descSC, SC, 1 );
        free( descSR );
        free( descSC );
    }

    CHAMELEON_zcesca_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zcesca_Tile_Async - Non-blocking equivalent of CHAMELEON_zcesca_Tile().
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
 * @sa CHAMELEON_zcesca
 * @sa CHAMELEON_zcesca_Tile
 * @sa CHAMELEON_scesca_Tile_Async
 *
 */
int CHAMELEON_zcesca_Tile_Async( int center, int scale, cham_store_t axis, CHAM_desc_t *A, void *user_ws,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    struct chameleon_pzcesca_s *ws;

    /* Check input arguments */
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zcesca_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( (center != 0) && (center != 1) ) {
        chameleon_error("CHAMELEON_zcesca_Tile_Async", "Illegal value of center");
        return -1;
    }
    if ( (scale != 0) && (scale != 1) ) {
        chameleon_error("CHAMELEON_zcesca_Tile_Async", "Illegal value of scale");
        return -2;
    }
    if ( (axis != ChamColumnwise) && (axis != ChamRowwise) && (axis != ChamEltwise) ) {
        chameleon_error("CHAMELEON_zcesca_Tile_Async", "Illegal value of axis");
        return -3;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zcesca_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zcesca_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zcesca_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zcesca_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zcesca_WS_Alloc( A );
    }
    else {
        ws = user_ws;
    }

    chameleon_pzcesca( ws, center, scale, axis, A, sequence, request );

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zcesca_WS_Free( ws );
    }

    return CHAMELEON_SUCCESS;
}
