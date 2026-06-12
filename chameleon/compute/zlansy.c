/**
 *
 * @file zlansy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
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
 *  CHAMELEON_zlansy returns the value
 *
 *     zlansy = ( max(abs(A(i,j))), NORM = ChamMaxNorm
 *              (
 *              ( norm1(A),         NORM = ChamOneNorm
 *              (
 *              ( normI(A),         NORM = ChamInfNorm
 *              (
 *              ( normF(A),         NORM = ChamFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The number of columns/rows of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval the norm described above.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlansy_Tile
 * @sa CHAMELEON_zlansy_Tile_Async
 * @sa CHAMELEON_clansy
 * @sa CHAMELEON_dlansy
 * @sa CHAMELEON_slansy
 *
 */
double CHAMELEON_zlansy(cham_normtype_t norm, cham_uplo_t uplo, int N,
                    CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    double value = -1.;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlansy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (norm != ChamMaxNorm) && (norm != ChamOneNorm)
         && (norm != ChamInfNorm) && (norm != ChamFrobeniusNorm) ) {
        chameleon_error("CHAMELEON_zlansy", "illegal value of norm");
        return -1;
    }
    if ( (uplo != ChamUpper) && (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlansy", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlansy", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zlansy", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if ( N == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlansy", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlansy_Tile_Async( norm, uplo, &descAt, &value, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return value;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zlansy_Tile - Tile equivalent of CHAMELEON_zlansy().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlansy
 * @sa CHAMELEON_zlansy_Tile_Async
 * @sa CHAMELEON_clansy_Tile
 * @sa CHAMELEON_dlansy_Tile
 * @sa CHAMELEON_slansy_Tile
 *
 */
double CHAMELEON_zlansy_Tile( cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;
    double value = -1.;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlansy_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlansy_Tile_Async( norm, uplo, A, &value, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return ( status == CHAMELEON_SUCCESS ) ? value : (double)status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zlansy_Tile_Async - Non-blocking equivalent of CHAMELEON_zlansy_Tile().
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
 * @sa CHAMELEON_zlansy
 * @sa CHAMELEON_zlansy_Tile
 * @sa CHAMELEON_clansy_Tile_Async
 * @sa CHAMELEON_dlansy_Tile_Async
 * @sa CHAMELEON_slansy_Tile_Async
 *
 */
int CHAMELEON_zlansy_Tile_Async( cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *value,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlansy_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlansy_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlansy_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zlansy_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zlansy_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (norm != ChamMaxNorm) && (norm != ChamOneNorm)
         && (norm != ChamInfNorm) && (norm != ChamFrobeniusNorm) ) {
        chameleon_error("CHAMELEON_zlansy_Tile", "illegal value of norm");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (uplo != ChamUpper) && (uplo != ChamLower) ) {
        chameleon_error("CHAMELEON_zlansy_Tile", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if ( A->m == 0) {
        *value = 0.0;
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlansy_generic( norm, uplo, ChamTrans, A, value, sequence, request );

    return CHAMELEON_SUCCESS;
}
