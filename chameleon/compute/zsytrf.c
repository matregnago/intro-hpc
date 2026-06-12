/**
 *
 * @file zsytrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf wrappers
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Marc Sergent
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zsytrf - Computes the Cholesky factorization of a symmetric matrix A.
 *  The factorization has the form
 *
 *    \f[ A = \{_{L\times L^T, if uplo = ChamLower}^{U^T\times U, if uplo = ChamUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U^H*U or A = L*L^H.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsytrf_Tile
 * @sa CHAMELEON_zsytrf_Tile_Async
 * @sa CHAMELEON_csytrf
 * @sa CHAMELEON_dsytrf
 * @sa CHAMELEON_ssytrf
 * @sa CHAMELEON_zpotrs
 *
 */
int CHAMELEON_zsytrf( cham_uplo_t uplo, int N,
                  CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsytrf", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zsytrf", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zsytrf", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zsytrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_max(N, 0) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZSYSV, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsytrf", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zsytrf_Tile_Async( uplo, &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
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
 *  CHAMELEON_zsytrf_Tile - Computes the Cholesky factorization of a symmetric matrix.
 *  Tile equivalent of CHAMELEON_zsytrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U^T*U or A = L*L^T.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsytrf
 * @sa CHAMELEON_zsytrf_Tile_Async
 * @sa CHAMELEON_csytrf_Tile
 *
 */
int CHAMELEON_zsytrf_Tile( cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsytrf_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zsytrf_Tile_Async( uplo, A, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zsytrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  matrix.
 *  Non-blocking equivalent of CHAMELEON_zsytrf_Tile().
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
 * @sa CHAMELEON_zsytrf
 * @sa CHAMELEON_zsytrf_Tile
 * @sa CHAMELEON_csytrf_Tile_Async
 *
 */
int CHAMELEON_zsytrf_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsytrf_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zsytrf_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zsytrf_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zsytrf_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zsytrf_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zsytrf_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -1);
    }
    /* Quick return */
    /*
     if (chameleon_max(N, 0) == 0)
     return CHAMELEON_SUCCESS;
     */
    chameleon_pzsytrf( uplo, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
