/**
 *
 * @file zlascal.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal wrappers
 *
 * @version 1.4.0
 * @author Dalal Sukkari
 * @author Mathieu Faverge
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
 *  CHAMELEON_zlascal - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = ChamUpperLower: A is a general matrix.
 *          = ChamUpper: A is an upper trapezoidal matrix.
 *          = ChamLower: A is a lower trapezoidal matrix.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in,out] A
 *          A is a LDA-by-N matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlascal_Tile
 * @sa CHAMELEON_clascal
 * @sa CHAMELEON_dlascal
 * @sa CHAMELEON_slascal
 *
 */
int CHAMELEON_zlascal( cham_uplo_t uplo, int M, int N,
                   CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    CHAM_desc_t descAl, descAt;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlascal", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower) && (uplo != ChamUpperLower)) {
        chameleon_error("CHAMELEON_zlascal", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zlascal", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlascal", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlascal", "illegal value of LDA");
        return -6;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        (alpha == (CHAMELEON_Complex64_t)1.0))
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlascal", "chameleon_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlascal_Tile_Async( uplo, alpha, &descAt, sequence, &request );

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
 *  CHAMELEON_zlascal_Tile - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = ChamUpperLower: A is a general matrix.
 *          = ChamUpper: A is an upper trapezoidal matrix.
 *          = ChamLower: A is a lower trapezoidal matrix.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-N matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlascal
 * @sa CHAMELEON_zlascal_Tile_Async
 * @sa CHAMELEON_clascal_Tile
 * @sa CHAMELEON_dlascal_Tile
 * @sa CHAMELEON_slascal_Tile
 *
 */
int CHAMELEON_zlascal_Tile( cham_uplo_t uplo,
                            CHAMELEON_Complex64_t alpha, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlascal_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlascal_Tile_Async( uplo, alpha, A, sequence, &request );

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
 *  CHAMELEON_zlascal_Tile_Async - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *  Non-blocking equivalent of CHAMELEON_zlascal_Tile().
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
 * @sa CHAMELEON_zlascal
 * @sa CHAMELEON_zlascal_Tile
 * @sa CHAMELEON_clascal_Tile_Async
 * @sa CHAMELEON_dlascal_Tile_Async
 * @sa CHAMELEON_slascal_Tile_Async
 *
 */
int CHAMELEON_zlascal_Tile_Async( cham_uplo_t uplo,
                              CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t descA;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlascal_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlascal_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlascal_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zlascal_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower) && (uplo != ChamUpperLower)) {
        chameleon_error("CHAMELEON_zlascal", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    if ( (descA.i%descA.mb != 0) || (descA.j%descA.nb != 0) ) {
        chameleon_error("CHAMELEON_zlascal", "start indexes have to be multiple of tile size");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (descA.m == 0) || (descA.n == 0) ||
         (alpha == (CHAMELEON_Complex64_t)1.0) )
    {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlascal( uplo, alpha, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
