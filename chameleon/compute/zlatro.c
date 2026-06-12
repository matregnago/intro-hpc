/**
 *
 * @file zlatro.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro wrappers
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2026-01-23
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zlatro transposes all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = ChamUpper, only the upper trapezium
 *          is accessed; if UPLO = ChamLower, only the lower trapezium is
 *          accessed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] B
 *          The N-by-M matrix B.
 *          On exit, B = A^T in the locations specified by UPLO.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatro_Tile
 * @sa CHAMELEON_zlatro_Tile_Async
 * @sa CHAMELEON_clatro
 * @sa CHAMELEON_dlatro
 * @sa CHAMELEON_slatro
 *
 */
int CHAMELEON_zlatro( cham_uplo_t uplo, cham_trans_t trans,
                      int M, int N,
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
    cham_uplo_t uploB;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatro", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlatro", "illegal value of uplo");
        return -1;
    }
    if( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "illegal value of trans");
        return -1;
    }
    /* quick return */
    if ( trans == ChamNoTrans ) {
        return CHAMELEON_SUCCESS;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zlatro", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlatro", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlatro", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zlatro", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlatro", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;
    uploB = uplo == ChamUpperLower ? ChamUpperLower :
                        ( uplo == ChamUpper ? ChamLower : ChamUpper );

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, uplo,
                         A, NB, NB, LDA, N, M, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, uploB,
                         B, NB, NB, LDB, M, N, M, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlatro_Tile_Async( uplo, trans, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, uplo, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                         ChamDescInout, uploB, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zlatro_Tile - Tile equivalent of CHAMELEON_zlatro().
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
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   A is not transposed;
 *          = ChamTrans:     A is transposed;
 *          = ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = ChamUpper, only the upper trapezium
 *          is accessed; if UPLO = ChamLower, only the lower trapezium is
 *          accessed.
 *
 * @param[out] B
 *          The N-by-M matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatro
 * @sa CHAMELEON_zlatro_Tile_Async
 * @sa CHAMELEON_clatro_Tile
 * @sa CHAMELEON_dlatro_Tile
 * @sa CHAMELEON_slatro_Tile
 *
 */
int CHAMELEON_zlatro_Tile( cham_uplo_t uplo, cham_trans_t trans,
                           CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatro_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlatro_Tile_Async( uplo, trans, A, B, sequence, &request );

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
 *  CHAMELEON_zlatro_Tile_Async - Non-blocking equivalent of CHAMELEON_zlatro_Tile().
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
 * @sa CHAMELEON_zlatro
 * @sa CHAMELEON_zlatro_Tile
 * @sa CHAMELEON_clatro_Tile_Async
 * @sa CHAMELEON_dlatro_Tile_Async
 * @sa CHAMELEON_slatro_Tile_Async
 *
 */
int CHAMELEON_zlatro_Tile_Async( cham_uplo_t uplo, cham_trans_t trans,
                                 CHAM_desc_t *A, CHAM_desc_t *B,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatro_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatro_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatro_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((A->mb != B->nb) || (A->nb != B->mb) ){
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "only matching tile sizes supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "illegal value of uplo");
        return -1;
    }
    if( !isValidTrans( trans ) ) {
        chameleon_error("CHAMELEON_zlatro_Tile_Async", "illegal value of trans");
        return -1;
    }
    /* quick return */
    if ( trans == ChamNoTrans ) {
        return CHAMELEON_SUCCESS;
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlatro( uplo, trans, A, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
