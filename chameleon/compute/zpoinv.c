/**
 *
 * @file zpoinv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpoinv wrappers
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Eyraud-Dubois Lionel
 * @author Florent Pruvost
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
 *  @brief Computes the inverse of a complex Hermitian positive definite
 *  matrix A using the Cholesky factorization A = U^H*U or A = L*L^H
 *  computed by CHAMELEON_zpotrf().
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
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the upper or lower triangle of the (Hermitian)
 *          inverse of A, overwriting the input factor U or L.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, the (i,i) element of the factor U or L is
 *                zero, and the inverse could not be computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zpoinv_Tile
 * @sa CHAMELEON_zpoinv_Tile_Async
 * @sa CHAMELEON_cpoinv
 * @sa CHAMELEON_dpoinv
 * @sa CHAMELEON_spoinv
 * @sa CHAMELEON_zpotrf
 *
 */
int CHAMELEON_zpoinv( cham_uplo_t uplo, int N,
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
        chameleon_fatal_error("CHAMELEON_zpoinv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zpoinv", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zpoinv", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zpoinv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_max(N, 0) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZPOSV, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zpoinv", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zpoinv_Tile_Async( uplo, &descAt, sequence, &request );

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
 * @brief Tile equivalent of CHAMELEON_zpoinv()
 *
 *  Computes the inverse of a complex Hermitian
 *  positive definite matrix A using the Cholesky factorization
 *  A = U^H*U or A = L*L^H computed by CHAMELEON_zpotrf().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  Can be further configured via the CHAMELEON_POINV_REMAP environment variable
 *  which specifies a custom data distribution file to be used for the TRTRI operation.
 *  In that case the order of operations is:
 *  1. POTRF
 *  2. redistribution to the provided custom distribution
 *  3. TRTRI
 *  4. redistribution back to the original distribution
 *  5. LAUUM
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the upper or lower triangle of the (Hermitian)
 *          inverse of A, overwriting the input factor U or L.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, the leading minor of order i of A is not
 *               positive definite, so the factorization could not be
 *               completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zpoinv
 * @sa CHAMELEON_zpoinv_Tile_Async
 * @sa CHAMELEON_cpoinv_Tile
 * @sa CHAMELEON_dpoinv_Tile
 * @sa CHAMELEON_spoinv_Tile
 * @sa CHAMELEON_zpotrf_Tile
 *
 */
int CHAMELEON_zpoinv_Tile( cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpoinv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zpoinv_Tile_Async( uplo, A, sequence, &request );

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
 * @brief Non-blocking equivalent of CHAMELEON_zpoinv_Tile().
 *
 * Submit the tasks to compute the inverse of a complex Hermitian
 * positive definite matrix A using the Cholesky factorization A = U^H*U
 * or A = L*L^H computed by CHAMELEON_zpotrf.
 *
 * May return before the computation is finished.
 * Allows for pipelining of operations at runtime.
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
 * @sa CHAMELEON_zpoinv
 * @sa CHAMELEON_zpoinv_Tile
 * @sa CHAMELEON_cpoinv_Tile_Async
 * @sa CHAMELEON_dpoinv_Tile_Async
 * @sa CHAMELEON_spoinv_Tile_Async
 * @sa CHAMELEON_zpotrf_Tile_Async
 *
 */
int CHAMELEON_zpoinv_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t  *chamctxt;
    int              flush, change_distribution_for_trtri = 0;
    custom_dist_t   *trtri_custom_get_rankof_arg;
    custom_dist_t   *original_get_rankof_arg;
    blkrankof_fct_t  original_get_rankof;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpoinv_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zpoinv_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zpoinv_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zpoinv_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zpoinv_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zpoinv_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -1);
    }
    /* Quick return */
    if ( chameleon_max( A->n, 0 ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    /* If the CHAMELEON_POINV_REMAP variable is present, use its value
       to create a custom distribution for the TRTRI operation */
    {
        char *custom_dist = chameleon_getenv( "CHAMELEON_POINV_REMAP" );
        int rc;

        if ( custom_dist != NULL ) {
            rc = chameleon_getrankof_custom_init( &trtri_custom_get_rankof_arg, custom_dist );
            if ( rc != CHAMELEON_SUCCESS ) {
                chameleon_error( "CHAMELEON_zpoinv", "CHAMELEON_POINV_REMAP: could not open the custom distribution file, keeping the original distribution");
            } else {
                change_distribution_for_trtri = 1;
                original_get_rankof     = A->get_rankof_init;
                original_get_rankof_arg = A->get_rankof_init_arg;
            }
        }

        chameleon_cleanenv( custom_dist );
    }

    flush = request->flush;
    request->flush = CHAMELEON_FALSE;
    chameleon_pzpotrf( uplo, A, sequence, request );

    if ( change_distribution_for_trtri ) {
        CHAMELEON_Desc_Change_Distribution_Async( uplo, A, chameleon_getrankof_custom, trtri_custom_get_rankof_arg, sequence );
    }

    chameleon_pztrtri( uplo, ChamNonUnit, A, sequence, request );

    if ( change_distribution_for_trtri ) {
        CHAMELEON_Desc_Change_Distribution_Async( uplo, A, original_get_rankof, original_get_rankof_arg, sequence );
    }

    request->flush = flush;
    chameleon_pzlauum( uplo, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
