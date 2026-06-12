/**
 *
 * @file zgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv wrappers
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
 *  CHAMELEON_zgetrf_incpiv - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descL
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by CHAMELEON_zgetrs_incpiv to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_incpiv_Tile
 * @sa CHAMELEON_zgetrf_incpiv_Tile_Async
 * @sa CHAMELEON_cgetrf_incpiv
 * @sa CHAMELEON_dgetrf_incpiv
 * @sa CHAMELEON_sgetrf_incpiv
 * @sa CHAMELEON_zgetrs_incpiv
 *
 */
int CHAMELEON_zgetrf_incpiv( int M, int N,
                         CHAMELEON_Complex64_t *A, int LDA,
                         CHAM_desc_t *descL, int *IPIV )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_incpiv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgetrf_incpiv", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgetrf_incpiv", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgetrf_incpiv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGESV, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrf_incpiv", "chameleon_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB   = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgetrf_incpiv_Tile_Async( &descAt, descL, IPIV, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descL, sequence );

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
 *  CHAMELEON_zgetrf_incpiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of CHAMELEON_zgetrf_incpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by CHAMELEON_zgetrs_incpiv to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_incpiv
 * @sa CHAMELEON_zgetrf_incpiv_Tile_Async
 * @sa CHAMELEON_cgetrf_incpiv_Tile
 * @sa CHAMELEON_dgetrf_incpiv_Tile
 * @sa CHAMELEON_sgetrf_incpiv_Tile
 * @sa CHAMELEON_zgetrs_incpiv_Tile
 *
 */
int CHAMELEON_zgetrf_incpiv_Tile( CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_incpiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgetrf_incpiv_Tile_Async( A, L, IPIV, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( L, sequence );

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
 *  CHAMELEON_zgetrf_incpiv_Tile_Async - Computes the tile LU factorization of a matrix.
 *  Non-blocking equivalent of CHAMELEON_zgetrf_incpiv_Tile().
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
 * @sa CHAMELEON_zgetrf_incpiv
 * @sa CHAMELEON_zgetrf_incpiv_Tile
 * @sa CHAMELEON_cgetrf_incpiv_Tile_Async
 * @sa CHAMELEON_dgetrf_incpiv_Tile_Async
 * @sa CHAMELEON_sgetrf_incpiv_Tile_Async
 * @sa CHAMELEON_zgetrs_incpiv_Tile_Async
 *
 */
int CHAMELEON_zgetrf_incpiv_Tile_Async( CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV,
                                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t D, *Dptr = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_incpiv_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_incpiv_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrf_incpiv_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgetrf_incpiv_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(L) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrf_incpiv_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zgetrf_incpiv_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(M, N) == 0)
     return CHAMELEON_SUCCESS;
     */

#if defined(CHAMELEON_COPY_DIAG)
    {
        int n = chameleon_min( A->m, A->n ) * A->nb;
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, A->m, n );
        Dptr = &D;
    }
#endif

    chameleon_pzgetrf_incpiv( A, L, Dptr, IPIV, sequence, request );

    if (Dptr != NULL) {
        /* Flush data to be able to sync with runtimes such as PaRSEC */
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( L, sequence );
        CHAMELEON_Desc_Flush( Dptr, sequence );

        /* Wait and cleanup */
        chameleon_sequence_wait( chamctxt, sequence );
        chameleon_desc_destroy( Dptr );
    }
    (void)D;
    return CHAMELEON_SUCCESS;
}
