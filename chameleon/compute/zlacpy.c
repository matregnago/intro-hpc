/**
 *
 * @file zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy wrappers
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
 *  CHAMELEON_zlacpy copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
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
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlacpy_Tile
 * @sa CHAMELEON_zlacpy_Tile_Async
 * @sa CHAMELEON_clacpy
 * @sa CHAMELEON_dlacpy
 * @sa CHAMELEON_slacpy
 *
 */
int CHAMELEON_zlacpy( cham_uplo_t uplo, int M, int N,
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

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlacpy", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zlacpy", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlacpy", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlacpy", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlacpy", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlacpy", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, uplo,
                     B, NB, NB, LDB, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlacpy_Tile_Async( uplo, &descAt, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, uplo, sequence, &request );

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
 *  CHAMELEON_zlacpy_Tile - Tile equivalent of CHAMELEON_zlacpy().
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
 * @param[in] A
 *          The M-by-N matrix A. If uplo = ChamUpper, only the upper trapezium
 *          is accessed; if UPLO = ChamLower, only the lower trapezium is
 *          accessed.
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlacpy
 * @sa CHAMELEON_zlacpy_Tile_Async
 * @sa CHAMELEON_clacpy_Tile
 * @sa CHAMELEON_dlacpy_Tile
 * @sa CHAMELEON_slacpy_Tile
 *
 */
int CHAMELEON_zlacpy_Tile( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlacpy_Tile_Async( uplo, A, B, sequence, &request );

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
 *  CHAMELEON_zlacpy_Tile_Async - Non-blocking equivalent of CHAMELEON_zlacpy_Tile().
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
 * @sa CHAMELEON_zlacpy
 * @sa CHAMELEON_zlacpy_Tile
 * @sa CHAMELEON_clacpy_Tile_Async
 * @sa CHAMELEON_dlacpy_Tile_Async
 * @sa CHAMELEON_slacpy_Tile_Async
 *
 */
int CHAMELEON_zlacpy_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zlacpy_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlacpy_Tile_Async", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((A->mb != B->mb) || (A->nb != B->nb) ){
        chameleon_error("CHAMELEON_zlacpy_Tile_Async", "only matching tile sizes supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper) &&
         (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlacpy_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlacpy( uplo, A, B, sequence, request );

    return CHAMELEON_SUCCESS;
}
