/**
 *
 * @file zgepdf_qdwh.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qdwh wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Computes the polar decomposition of a general M-by-N matrix A = U_p H.
 *
 * \[ A = U_p * H \] where:
 *     - U_p is M-by-N and \[ U_p U_p^t = I \]
 *     - H is N-by-N and is hermitian/symmetric positive semidefinite.
 *
 * On output, A is overwritten by U.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrices A and H, and number of rows of
 *          the matrix H. N >= 0.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of A2.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in,out] A
 *          On entry, the LDA-by-N matrix to factorize.
 *          On exit, contains the orthonormal U_p matrix of the polar decomposition.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] H
 *          On exit, the LDH-by-N matrix that holds the polar factor of the
 *          decomposition of A.
 *
 * @param[in,out] info
 *          On entry, the allocated data structure to store informations on the
 *          polar decomposition.
 *          On exit, contains the number of each iterations: QR and Cholesky, as
 *          well as the number of flops of the algorithm.
 *          Not referenced, if null.
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgepdf_qdwh_Tile
 * @sa CHAMELEON_zgepdf_qdwh_Tile_Async
 * @sa CHAMELEON_cgepdf_qdwh
 * @sa CHAMELEON_dgepdf_qdwh
 * @sa CHAMELEON_sgepdf_qdwh
 *
 */
int CHAMELEON_zgepdf_qdwh( int M, int N,
                           CHAMELEON_Complex64_t *A, int LDA,
                           CHAMELEON_Complex64_t *H, int LDH,
                           gepdf_info_t *info )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descHl, descHt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qdwh", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgepdf_qdwh", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgepdf_qdwh", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgepdf_qdwh", "illegal value of LDA");
        return -4;
    }
    if (LDH < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zgepdf_qdwh", "illegal value of LDH");
        return -6;
    }

    /* Quick return */
    if ( (M == 0) || (N == 0) ) {
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgepdf_qdwh", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "H", &descHl, &descHt, ChamDescOutput, ChamUpperLower,
                         H, NB, NB, LDH, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgepdf_qdwh_Tile_Async( &descAt, &descHt, info,
                                      sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descHl, &descHt,
                         ChamDescOutput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descHl, &descHt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  @brief Compute the polar decomposition of A = U_p * H.
 *  Tile equivalent of CHAMELEON_zgepdf_qdwh().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry the M-by-N input matrix A for which the polar decomposition
 *          is computed.
 *          On exit, A is overwritten by U_p.
 *
 * @param[out] H
 *          Not referenced on entry.
 *          On exit, contains the symmetric/hermitian polar factor.
 *
 * @param[in,out] info
 *          On entry, the allocated data structure to store informations on the
 *          polar decomposition.
 *          On exit, contains the number of each iterations: QR and Cholesky, as
 *          well as the number of flops of the algorithm.
 *          Not referenced, if null.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS on successful exit
 * @retval <0 on failure
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgepdf_qdwh
 * @sa CHAMELEON_zgepdf_qdwh_Tile_Async
 * @sa CHAMELEON_cgepdf_qdwh_Tile
 * @sa CHAMELEON_dgepdf_qdwh_Tile
 * @sa CHAMELEON_sgepdf_qdwh_Tile
 *
 */
int CHAMELEON_zgepdf_qdwh_Tile( CHAM_desc_t *A, CHAM_desc_t *H, gepdf_info_t *info )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qdwh_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgepdf_qdwh_Tile_Async( A, H, info,
                                      sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( H, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  @brief Compute the polar decomposition of A = U_p * H.
 *  Non-blocking equivalent of CHAMELEON_zgepdf_qdwh_Tile().
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
 * @sa CHAMELEON_zgepdf_qdwh
 * @sa CHAMELEON_zgepdf_qdwh_Tile
 * @sa CHAMELEON_cgepdf_qdwh_Tile_Async
 * @sa CHAMELEON_dgepdf_qdwh_Tile_Async
 * @sa CHAMELEON_sgepdf_qdwh_Tile_Async
 *
 */
int CHAMELEON_zgepdf_qdwh_Tile_Async( CHAM_desc_t *A, CHAM_desc_t *H, gepdf_info_t *info,
                                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_zgepdf_qdwh_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qdwh_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qdwh_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    if ( (A->mb != A->nb) || (H->mb != H->nb) || (A->mb != H->mb) ) {
        chameleon_error("CHAMELEON_zgepdf_qdwh_Tile_Async", "tile sizes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (A->n != H->m) || (A->n != H->n) ) {
        chameleon_error("CHAMELEON_zgepdf_qdwh_Tile_Async", "H must be square and of size the number of column of A");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (A->m == 0) || (A->n == 0) ) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzgepdf_qdwh( ChamGeneral, A, H,
                            info, sequence, request );

    return CHAMELEON_SUCCESS;
}
