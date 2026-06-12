/**
 *
 * @file ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm wrappers
 *
 * @version 1.4.0
 * @author Jakub Kurzak
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
 *  CHAMELEON_ztrsm - Solves one of the matrix equations
 *    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
 *
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
 *
 *    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
 *
 * The matrix X is overwritten on B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether op(A) appears on the left or on the right of X:
 *          = ChamLeft:  op(A) * X = B
 *          = ChamRight: X * op(A) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: A is an upper triangular matrix.
 *          = ChamLower: A is a lower triangular matrix.
 *
 * @param[in] trans
 *          Specifies the form of op( A ) to be used in the matrix
 *          multiplication as follows:
 *          = ChamNoTrans:   op( A ) = A.
 *          = ChamTrans:     op( A ) = A**T.
 *          = ChamConjTrans: op( A ) = A**H.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: A is not assumed to be unit triangular.
 *          = ChamUnit:    A is assumed to be unit triangular.
 *
 * @param[in] M
 *          M specifies the number of rows of B. M must be at least zero.
 *
 * @param[in] N
 *          N specifies the number of columns of B. N must be at least zero.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha. When alpha is zero then A is not
 *          referenced and B need not be set before entry.
 *
 * @param[in] A
 *          The triangular matrix A of dimension lda-by-k,
 *          where k is M when side = ChamLeft
 *            and k is N when side = ChamRight
 *          Before entry with uplo = ChamUpper, the leading k-by-k upper
 *          triangular part of the array A must contain the upper triangular
 *          matrix and the strictly lower triangular part of A is not
 *          referenced.
 *          Before entry with uplo = ChamLower, the leading k-by-k lower
 *          triangular part of the array A must contain the lower triangular
 *          matrix and the strictly upper triangular part of A is not
 *          referenced.
 *          Note that when diag = ChamUnit, the diagonal elements of A are not
 *          referenced either, but are assumed to be unity.
 *
 * @param[in] LDA
 *          LDA specifies the first dimension of A. When side = ChamLeft then
 *          LDA must be at least max( 1, M ), when side = ChamRight then LDA
 *          must be at least max( 1, N ).
 *
 * @param[in,out] B
 *          The matrix B of dimension LDB-by-N.
 *          Before entry, the leading m by n part of the array B must contain
 *          the right-hand side matrix B, and on exit is overwritten by the
 *          solution matrix X.
 *
 * @param[in] LDB
 *          LDB specifies the first dimension of B. LDB must be at least
 *          max( 1, M ).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrsm_Tile
 * @sa CHAMELEON_ztrsm_Tile_Async
 * @sa CHAMELEON_ctrsm
 * @sa CHAMELEON_dtrsm
 * @sa CHAMELEON_strsm
 *
 */
int CHAMELEON_ztrsm( cham_side_t side, cham_uplo_t uplo,
                     cham_trans_t trans, cham_diag_t diag,
                     int M, int N, CHAMELEON_Complex64_t alpha,
                     CHAMELEON_Complex64_t *A, int LDA,
                     CHAMELEON_Complex64_t *B, int LDB )
{
    int NB, Ak;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ( side == ChamLeft ) {
        Ak = M;
    } else {
        Ak = N;
    }

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of side");
        return -1;
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of uplo");
        return -2;
    }
    if ( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of trans");
        return -3;
    }
    if ((diag != ChamUnit) && (diag != ChamNonUnit)) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of diag");
        return -4;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of M");
        return -5;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of N");
        return -6;
    }
    if (LDA < chameleon_max(1, Ak)) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of LDA");
        return -8;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_ztrsm", "illegal value of LDB");
        return -10;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, M & N; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZPOSV, M, M, N);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztrsm", "chameleon_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, Ak, Ak, Ak, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, N,  M,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_ztrsm_Tile_Async(  side, uplo, trans, diag, alpha, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );
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
 *  CHAMELEON_ztrsm_Tile - Solves one of the matrix equations
 *    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B.
 *
 *  Tile equivalent of CHAMELEON_ztrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether op(A) appears on the left or on the right of X:
 *          = ChamLeft:  op(A) * X = B
 *          = ChamRight: X * op(A) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: A is an upper triangular matrix.
 *          = ChamLower: A is a lower triangular matrix.
 *
 * @param[in] trans
 *          Specifies the form of op( A ) to be used in the matrix
 *          multiplication as follows:
 *          = ChamNoTrans:   op( A ) = A.
 *          = ChamTrans:     op( A ) = A**T.
 *          = ChamConjTrans: op( A ) = A**H.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: A is assumed to be unit triangular.
 *          = ChamUnit:    A is not assumed to be unit triangular.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha. When alpha is zero then A is not
 *          referenced and B need not be set before entry.
 *
 * @param[in] A
 *          The triangular matrix A of dimension lda-by-k,
 *          where k is M when side = ChamLeft
 *            and k is N when side = ChamRight
 *          Before entry with uplo = ChamUpper, the leading k-by-k upper
 *          triangular part of the array A must contain the upper triangular
 *          matrix and the strictly lower triangular part of A is not
 *          referenced.
 *          Before entry with uplo = ChamLower, the leading k-by-k lower
 *          triangular part of the array A must contain the lower triangular
 *          matrix and the strictly upper triangular part of A is not
 *          referenced.
 *          Note that when diag = ChamUnit, the diagonal elements of A are not
 *          referenced either, but are assumed to be unity.
 *
 * @param[in,out] B
 *          The matrix B of dimension LDB-by-N.
 *          Before entry, the leading m by n part of the array B must contain
 *          the right-hand side matrix B, and on exit is overwritten by the
 *          solution matrix X.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrsm
 * @sa CHAMELEON_ztrsm_Tile_Async
 * @sa CHAMELEON_ctrsm_Tile
 * @sa CHAMELEON_dtrsm_Tile
 * @sa CHAMELEON_strsm_Tile
 *
 */
int CHAMELEON_ztrsm_Tile( cham_side_t side, cham_uplo_t uplo,
                          cham_trans_t trans, cham_diag_t diag,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_ztrsm_Tile_Async(side, uplo, trans, diag, alpha, A, B, sequence, &request );

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
 *  CHAMELEON_ztrsm_Tile_Async - Computes triangular solve.
 *  Non-blocking equivalent of CHAMELEON_ztrsm_Tile().
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
 * @sa CHAMELEON_ztrsm
 * @sa CHAMELEON_ztrsm_Tile
 * @sa CHAMELEON_ctrsm_Tile_Async
 * @sa CHAMELEON_dtrsm_Tile_Async
 * @sa CHAMELEON_strsm_Tile_Async
 *
 */
int CHAMELEON_ztrsm_Tile_Async( cham_side_t side, cham_uplo_t uplo,
                                cham_trans_t trans, cham_diag_t diag,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_ztrsm_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((A->nb != A->mb) || (B->nb != B->mb)) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (side != ChamLeft && side != ChamRight) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "illegal value of side");
        return chameleon_request_fail(sequence, request, -1);
    }
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -2);
    }
    if ( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "illegal value of trans");
        return chameleon_request_fail(sequence, request, -3);
    }
    if ((diag != ChamUnit) && (diag != ChamNonUnit)) {
        chameleon_error("CHAMELEON_ztrsm_Tile", "illegal value of diag");
        return chameleon_request_fail(sequence, request, -4);
    }

    /* Quick return */
    chameleon_pztrsm( side, uplo, trans, diag, alpha,
                      A, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
