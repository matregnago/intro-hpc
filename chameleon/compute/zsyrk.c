/**
 *
 * @file zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk wrappers
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
 *  CHAMELEON_zsyrk - Performs one of the hermitian rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of C is stored;
 *          = ChamLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans  :   A is transposed.
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = ChamNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ) if trans == ChamNoTrans, otherwise LDA must
 *          be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsyrk_Tile
 * @sa CHAMELEON_csyrk
 * @sa CHAMELEON_dsyrk
 * @sa CHAMELEON_ssyrk
 *
 */
int CHAMELEON_zsyrk( cham_uplo_t uplo, cham_trans_t trans, int N, int K,
                 CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t beta,  CHAMELEON_Complex64_t *C, int LDC )
{
    int NB;
    int Am, An;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descCl, descCt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsyrk", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of uplo");
        return -1;
    }
    if ((trans != ChamNoTrans) && (trans != ChamTrans)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of trans");
        return -2;
    }
    if ( trans == ChamNoTrans ) {
        Am = N; An = K;
    } else {
        Am = K; An = N;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of N");
        return -3;
    }
    if (K < 0) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of K");
        return -4;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of LDA");
        return -7;
    }
    if (LDC < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of LDC");
        return -10;
    }

    /* Quick return */
    if (N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0 || K == 0.0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZSYRK, N, K, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsyrk", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                     A, NB, NB, LDA, An, Am, An, sequence, &request );
    chameleon_zlap2tile( chamctxt, "C", &descCl, &descCt, ChamDescInout, uplo,
                     C, NB, NB, LDC, N, N,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zsyrk_Tile_Async( uplo, trans, alpha, &descAt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                     ChamDescInout, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
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
 *  CHAMELEON_zsyrk_Tile - Performs rank k update.
 *  Tile equivalent of CHAMELEON_zsyrk().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of C is stored;
 *          = ChamLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = ChamNoTrans,
 *          and is N otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zsyrk_Tile
 * @sa CHAMELEON_csyrk
 * @sa CHAMELEON_dsyrk
 * @sa CHAMELEON_ssyrk
 *
 */
int CHAMELEON_zsyrk_Tile( cham_uplo_t uplo, cham_trans_t trans,
                      CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                      CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsyrk_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zsyrk_Tile_Async( uplo, trans, alpha, A, beta, C, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
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
 *  CHAMELEON_zsyrk_Tile_Async - Performs rank-k update.
 *  Non-blocking equivalent of CHAMELEON_zsyrk_Tile().
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
 * @sa CHAMELEON_zsyrk
 * @sa CHAMELEON_zsyrk_Tile
 * @sa CHAMELEON_csyrk_Tile_Async
 * @sa CHAMELEON_dsyrk_Tile_Async
 * @sa CHAMELEON_ssyrk_Tile_Async
 *
 */
int CHAMELEON_zsyrk_Tile_Async( cham_uplo_t uplo, cham_trans_t trans,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                            CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    int N, K;
    int Am, An, Amb;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zsyrk_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zsyrk_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zsyrk_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -1);
    }
    if ((trans != ChamNoTrans) && (trans != ChamTrans)) {
        chameleon_error("CHAMELEON_zsyrk", "illegal value of transA");
        return chameleon_request_fail(sequence, request, -2);
    }

    if ( trans == ChamNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
    }

    if (C->mb != C->nb) {
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (Amb != C->mb) {
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "tile sizes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (C->m != C->n) {
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "only square matrix C is supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (Am != C->m) {
        chameleon_error("CHAMELEON_zsyrk_Tile_Async", "sizes of matrices have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    N = C->m;
    K = An;

    /* Quick return */
    if ( N == 0 ||
         ((alpha == (CHAMELEON_Complex64_t)0.0 || K == 0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    chameleon_pzsyrk( uplo, trans, alpha, A, beta, C, sequence, request );

    return CHAMELEON_SUCCESS;
}
