/**
 *
 * @file zgenm2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgenm2 wrappers
 *
 * @version 1.4.0
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
 * @brief Returns an estimate of the two-norm of A.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          Specify the tolerance used in the iterative algorithm to converge to
 *          an estimation. Must be > 0.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0. When M = 0,
 *          the returned value is set to zero.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval the two-norm estimate.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgenm2_Tile
 * @sa CHAMELEON_zgenm2_Tile_Async
 * @sa CHAMELEON_clange
 * @sa CHAMELEON_dlange
 * @sa CHAMELEON_slange
 *
 */
double CHAMELEON_zgenm2( double tol, int M, int N,
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
        chameleon_fatal_error("CHAMELEON_zgenm2", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgenm2", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgenm2", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgenm2", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgenm2", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgenm2_Tile_Async( tol, &descAt, &value, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );

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
 * @brief Tile equivalent of CHAMELEON_zgenm2().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          Specify the tolerance used in the iterative algorithm to converge to
 *          an estimation. Must be > 0.
 *
 * @param[in] A
 *          On entry, the input matrix A.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgenm2
 * @sa CHAMELEON_zgenm2_Tile_Async
 * @sa CHAMELEON_clange_Tile
 * @sa CHAMELEON_dlange_Tile
 * @sa CHAMELEON_slange_Tile
 *
 */
double CHAMELEON_zgenm2_Tile( double tol, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;
    double value = -1.;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgenm2_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgenm2_Tile_Async( tol, A, &value, sequence, &request );

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
 * @brief Non-blocking equivalent of CHAMELEON_zgenm2_Tile().
 *
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
 * @sa CHAMELEON_zgenm2
 * @sa CHAMELEON_zgenm2_Tile
 * @sa CHAMELEON_clange_Tile_Async
 * @sa CHAMELEON_dlange_Tile_Async
 * @sa CHAMELEON_slange_Tile_Async
 *
 */
int CHAMELEON_zgenm2_Tile_Async( double tol, CHAM_desc_t *A, double *value,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgenm2_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgenm2_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgenm2_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgenm2_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zgenm2_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( tol <= 0. ) {
        chameleon_error("CHAMELEON_zgenm2_Tile_Async", "tol must be > 0");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        *value = 0.0;
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzgenm2( tol, A, value, sequence, request );

    return CHAMELEON_SUCCESS;
}
