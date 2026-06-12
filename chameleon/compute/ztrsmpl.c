/**
 *
 * @file ztrsmpl.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsmpl wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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
 *  CHAMELEON_ztrsmpl - Performs the forward substitution step of solving a system of linear equations
 *  after the tile LU factorization of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by CHAMELEON_zgetrf_incpiv.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] descL
 *          Auxiliary factorization data, related to the tile L factor, computed by CHAMELEON_zgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from CHAMELEON_zgetrf_incpiv (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrsmpl_Tile
 * @sa CHAMELEON_ztrsmpl_Tile_Async
 * @sa CHAMELEON_ctrsmpl
 * @sa CHAMELEON_dtrsmpl
 * @sa CHAMELEON_strsmpl
 * @sa CHAMELEON_zgetrf_incpiv
 *
 */
int CHAMELEON_ztrsmpl( int N, int NRHS,
                   CHAMELEON_Complex64_t *A, int LDA,
                   CHAM_desc_t *descL, int *IPIV,
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
        chameleon_fatal_error("CHAMELEON_ztrsmpl", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        chameleon_error("CHAMELEON_ztrsmpl", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        chameleon_error("CHAMELEON_ztrsmpl", "illegal value of NRHS");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_ztrsmpl", "illegal value of LDA");
        return -4;
    }
    if (LDB < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_ztrsmpl", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(N, NRHS) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGESV, N, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztrsmpl", "chameleon_tune() failed");
        return status;
    }

    /* Set Mt, NT & NTRHS */
    NB    = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamLower,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_ztrsmpl_Tile_Async( &descAt, descL, IPIV, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descL, sequence );

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
 * CHAMELEON_ztrsmpl_Tile - Performs the forward substitution step of solving a system of linear equations
 * after the tile LU factorization of the matrix.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factor L from the factorization, computed by CHAMELEON_zgetrf_incpiv.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by CHAMELEON_zgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from CHAMELEON_zgetrf_incpiv (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztrsmpl
 * @sa CHAMELEON_ztrsmpl_Tile_Async
 * @sa CHAMELEON_ctrsmpl_Tile
 * @sa CHAMELEON_dtrsmpl_Tile
 * @sa CHAMELEON_strsmpl_Tile
 * @sa CHAMELEON_zgetrf_incpiv_Tile
 *
 */
int CHAMELEON_ztrsmpl_Tile( CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsmpl_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_ztrsmpl_Tile_Async( A, L, IPIV, B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( L, sequence );
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
 *  CHAMELEON_ztrsmpl_Tile - Performs the forward substitution step of solving
 *  a system of linear equations after the tile LU factorization of the matrix.
 *  Non-blocking equivalent of CHAMELEON_ztrsmpl_Tile().
 *  Returns control to the user thread before worker threads finish the computation
 *  to allow for pipelined execution of diferent routines.
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
 * @sa CHAMELEON_ztrsmpl
 * @sa CHAMELEON_ztrsmpl_Tile
 * @sa CHAMELEON_ctrsmpl_Tile_Async
 * @sa CHAMELEON_dtrsmpl_Tile_Async
 * @sa CHAMELEON_strsmpl_Tile_Async
 * @sa CHAMELEON_zgetrf_incpiv_Tile_Async
 *
 */
int CHAMELEON_ztrsmpl_Tile_Async( CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsmpl_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsmpl_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsmpl_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_ztrsmpl_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(L) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztrsmpl_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztrsmpl_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        chameleon_error("CHAMELEON_ztrsmpl_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(N, NRHS) == 0)
     return CHAMELEON_SUCCESS;
     */
    chameleon_pztrsmpl( A, B, L, IPIV, sequence, request );

    return CHAMELEON_SUCCESS;
}
