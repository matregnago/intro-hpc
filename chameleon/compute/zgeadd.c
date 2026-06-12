/**
 *
 * @file zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
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
 *  CHAMELEON_zgeadd - Performs a matrix addition similarly to the pzgeadd()
 *  function from the PBLAS library:
 *
 *    \f[ B = \alpha op( A ) + \beta B \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, and B are matrices, with op( A ) and B
 *  two m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix B. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( A ) and of the matrix B. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is N when trans = ChamNoTrans,
 *          and is M otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K), where K is M
 *          when trans = ChamNoTrans, and is N when otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDB-by-N matrix.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgeadd_Tile
 * @sa CHAMELEON_cgeadd
 * @sa CHAMELEON_dgeadd
 * @sa CHAMELEON_sgeadd
 *
 */
int CHAMELEON_zgeadd( cham_trans_t trans, int M, int N,
                      CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                      CHAMELEON_Complex64_t beta,  CHAMELEON_Complex64_t *B, int LDB )
{
    int NB;
    int Am, An;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgeadd", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_zgeadd", "illegal value of trans");
        return -1;
    }
    if ( trans == ChamNoTrans ) {
        Am = M; An = N;
    } else {
        Am = N; An = M;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zgeadd", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgeadd", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zgeadd", "illegal value of LDA");
        return -6;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgeadd", "illegal value of LDB");
        return -9;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (CHAMELEON_Complex64_t)0.0) && beta == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgeadd", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                     A, NB, NB, LDA, An, Am, An, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgeadd_Tile_Async( trans, alpha, &descAt, beta, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
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
 *  CHAMELEON_zgeadd_Tile - Performs a matrix addition similarly to the pzgeadd()
 *  function from the PBLAS library.
 *  Tile equivalent of CHAMELEON_zgeadd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is N when trans = ChamNoTrans,
 *          and is M otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDB-by-N matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgeadd
 * @sa CHAMELEON_zgeadd_Tile_Async
 * @sa CHAMELEON_cgeadd_Tile
 * @sa CHAMELEON_dgeadd_Tile
 * @sa CHAMELEON_sgeadd_Tile
 *
 */
int CHAMELEON_zgeadd_Tile( cham_trans_t trans,
                           CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                           CHAMELEON_Complex64_t beta,  CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgeadd_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgeadd_Tile_Async( trans, alpha, A, beta, B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

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
 *  CHAMELEON_zgeadd_Tile_Async - Performs a matrix addition similarly to the
 *  pzgeadd() function from the PBLAS library.
 *  Non-blocking equivalent of CHAMELEON_zgeadd_Tile().
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
 * @sa CHAMELEON_zgeadd
 * @sa CHAMELEON_zgeadd_Tile
 * @sa CHAMELEON_cgeadd_Tile_Async
 * @sa CHAMELEON_dgeadd_Tile_Async
 * @sa CHAMELEON_sgeadd_Tile_Async
 *
 */
int CHAMELEON_zgeadd_Tile_Async( cham_trans_t trans,
                                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                                 CHAMELEON_Complex64_t beta,  CHAM_desc_t *B,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    int M, N;
    int Am, An, Ai, Aj, Amb, Anb;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgeadd_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgeadd_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgeadd_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "illegal value of trans");
        return chameleon_request_fail(sequence, request, -1);
    }

    if ( trans == ChamNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }

    if ( (Amb != B->mb) || (Anb != B->nb) ) {
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "tile sizes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != B->m) || (An != B->n) ) {
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "sizes of matrices have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != B->i) || (Aj != B->j) ) {
        chameleon_error("CHAMELEON_zgeadd_Tile_Async", "start indexes have to match");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    M = B->m;
    N = B->n;

    /* Quick return */
    if ( (M == 0) || (N == 0) ||
         ((alpha == (CHAMELEON_Complex64_t)0.0) && (beta == (CHAMELEON_Complex64_t)1.0)) )
    {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pztradd( ChamUpperLower, trans, alpha, A, beta, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
