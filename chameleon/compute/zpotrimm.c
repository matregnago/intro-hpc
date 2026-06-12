/**
 *
 * @file zpotrimm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrimm wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
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
 *  CHAMELEON_zpotrimm - Computes the inverse of a complex Hermitian positive definite
 *  matrix A using the Cholesky factorization A = U^H*U or A = L*L^H
 *  computed by CHAMELEON_zpotrf.
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
 *          On entry, the triangular factor U or L from the Cholesky
 *          factorization A = U^H*U or A = L*L^H, as computed by
 *          CHAMELEON_zpotrf.
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
 * @sa CHAMELEON_zpotrimm_Tile
 * @sa CHAMELEON_zpotrimm_Tile_Async
 * @sa CHAMELEON_cpotri
 * @sa CHAMELEON_dpotri
 * @sa CHAMELEON_spotri
 * @sa CHAMELEON_zpotrf
 *
 */
int CHAMELEON_zpotrimm( cham_uplo_t uplo, int N,
                    CHAMELEON_Complex64_t *A, int LDA,
                    CHAMELEON_Complex64_t *B, int LDB,
                    CHAMELEON_Complex64_t *C, int LDC )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    CHAM_desc_t descCl, descCt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrimm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zpotrimm", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zpotrimm", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zpotrimm", "illegal value of LDA");
        return -4;
    }
    if (LDB < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zpotrimm", "illegal value of LDB");
        return -6;
    }
    if (LDC < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zpotrimm", "illegal value of LDC");
        return -8;
    }
    /* Quick return */
    if (chameleon_max(N, 0) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZPOSV, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zpotrimm", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInput, ChamUpperLower,
                     B, NB, NB, LDB, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "C", &descCl, &descCt, ChamDescInout, ChamUpperLower,
                     C, NB, NB, LDC, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zpotrimm_Tile_Async( uplo, &descAt, &descBt, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, uplo, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );
    chameleon_ztile2lap_cleanup( chamctxt, &descCl, &descCt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zpotrimm_Tile - Computes the inverse of a complex Hermitian
 *  positive definite matrix A using the Cholesky factorization
 *  A = U^H*U or A = L*L^H computed by CHAMELEON_zpotrf.
 *  Tile equivalent of CHAMELEON_zpotrimm().
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
 *          On entry, the triangular factor U or L from the Cholesky
 *          factorization A = U^H*U or A = L*L^H, as computed by
 *          CHAMELEON_zpotrf.
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
 * @sa CHAMELEON_zpotrimm
 * @sa CHAMELEON_zpotrimm_Tile_Async
 * @sa CHAMELEON_cpotri_Tile
 * @sa CHAMELEON_dpotri_Tile
 * @sa CHAMELEON_spotri_Tile
 * @sa CHAMELEON_zpotrf_Tile
 *
 */
int CHAMELEON_zpotrimm_Tile( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrimm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zpotrimm_Tile_Async( uplo, A, B, C, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

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
 *  CHAMELEON_zpotrimm_Tile_Async - Computes the inverse of a complex Hermitian
 *  positive definite matrix A using the Cholesky factorization A = U^H*U
 *  or A = L*L^H computed by CHAMELEON_zpotrf.
 *  Non-blocking equivalent of CHAMELEON_zpotrimm_Tile().
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
 * @sa CHAMELEON_zpotrimm
 * @sa CHAMELEON_zpotrimm_Tile
 * @sa CHAMELEON_cpotri_Tile_Async
 * @sa CHAMELEON_dpotri_Tile_Async
 * @sa CHAMELEON_spotri_Tile_Async
 * @sa CHAMELEON_zpotrf_Tile_Async
 *
 */
int CHAMELEON_zpotrimm_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrimm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrimm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrimm_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (B->nb != B->mb) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (C->nb != C->mb) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zpotrimm_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -1);
    }
    /* Quick return */
    /*
     if (chameleon_max(N, 0) == 0)
     return CHAMELEON_SUCCESS;
     */
    chameleon_pzpotrimm( uplo, A, B, C, sequence, request );
    /*
     chameleon_pztrtri( uplo, ChamNonUnit, A, sequence, request );
     chameleon_pzlauum( uplo, A, sequence, request );
     */

    return CHAMELEON_SUCCESS;
}
