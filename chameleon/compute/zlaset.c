/**
 *
 * @file zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset wrappers
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
 *  CHAMELEON_zlaset copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part is set. The lower
 *            triangle is unchanged.
 *            = ChamLower: Lower triangular part is set. The upper
 *            triangle is unchange.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          All the offdiagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaset_Tile
 * @sa CHAMELEON_zlaset_Tile_Async
 * @sa CHAMELEON_claset
 * @sa CHAMELEON_dlaset
 * @sa CHAMELEON_slaset
 *
 ******************************************************************************/
int CHAMELEON_zlaset( cham_uplo_t uplo, int M, int N,
                      CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
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
        chameleon_fatal_error("CHAMELEON_zlaset", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlaset", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zlaset", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlaset", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlaset", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlaset", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlaset_Tile_Async( uplo, alpha, beta, &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zlaset_Tile - Tile equivalent of CHAMELEON_zlaset().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaset
 * @sa CHAMELEON_zlaset_Tile_Async
 * @sa CHAMELEON_claset_Tile
 * @sa CHAMELEON_dlaset_Tile
 * @sa CHAMELEON_slaset_Tile
 *
 ******************************************************************************/
int CHAMELEON_zlaset_Tile( cham_uplo_t uplo,
                           CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                           CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlaset_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlaset_Tile_Async( uplo, alpha, beta, A, sequence, &request );

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
 *  CHAMELEON_zlaset_Tile_Async - Non-blocking equivalent of CHAMELEON_zlaset_Tile().
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
 * @sa CHAMELEON_zlaset
 * @sa CHAMELEON_zlaset_Tile
 * @sa CHAMELEON_claset_Tile_Async
 * @sa CHAMELEON_dlaset_Tile_Async
 * @sa CHAMELEON_slaset_Tile_Async
 *
 ******************************************************************************/
int CHAMELEON_zlaset_Tile_Async( cham_uplo_t uplo,
                                 CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                                 CHAM_desc_t *A,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlaset_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlaset_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlaset_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zlaset_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (alpha != beta) && (A->nb != A->mb) ) {
        chameleon_error("CHAMELEON_zlaset_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlaset_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlaset( uplo, alpha, beta, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
