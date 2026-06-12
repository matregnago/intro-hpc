/**
 *
 * @file zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
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
 * @brief Return the norm of A.
 *
 * CHAMELEON_zlange returns the value
 *     zlange = ( max(abs(A(i,j))), NORM = ChamMaxNorm
 *              (
 *              ( norm1(A),         NORM = ChamOneNorm
 *              (
 *              ( normI(A),         NORM = ChamInfNorm
 *              (
 *              ( normF(A),         NORM = ChamFrobeniusNorm
 *              (
 *              ( norm2(A),         NORM = ChamTwoNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum), normI
 *  denotes the infinity norm of a matrix (maximum row sum), normF denotes the
 *  Frobenius norm of a matrix (square root of sum of squares) and norm2 denotes
 *  an estimator of the maximum singular value. Note that max(abs(A(i,j))) is
 *  not a consistent matrix norm.
 *
 *  Note that ChamTwoNorm returns an estimate as computed by CHAMELEON_zgenm2()
 *  with a tolerance of 1.e-1.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
 *          = ChamTwoNorm: Twonorm
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
 * @retval the norm described above.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlange_Tile
 * @sa CHAMELEON_zlange_Tile_Async
 * @sa CHAMELEON_clange
 * @sa CHAMELEON_dlange
 * @sa CHAMELEON_slange
 *
 */
double CHAMELEON_zlange( cham_normtype_t norm, int M, int N,
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
        chameleon_fatal_error("CHAMELEON_zlange", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (norm != ChamMaxNorm) &&
         (norm != ChamOneNorm) &&
         (norm != ChamInfNorm) &&
         (norm != ChamFrobeniusNorm) &&
         (norm != ChamTwoNorm) )
    {
        chameleon_error("CHAMELEON_zlange", "illegal value of norm");
        return -1;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zlange", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlange", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlange", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0)
        return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlange", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlange_Tile_Async( norm, &descAt, &value, sequence, &request );

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
 * @brief Tile equivalent of CHAMELEON_zlange().
 *
 * Operates on matrices stored by tiles.
 * All matrices are passed through descriptors.
 * All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
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
 * @sa CHAMELEON_zlange
 * @sa CHAMELEON_zlange_Tile_Async
 * @sa CHAMELEON_clange_Tile
 * @sa CHAMELEON_dlange_Tile
 * @sa CHAMELEON_slange_Tile
 *
 */
double CHAMELEON_zlange_Tile( cham_normtype_t norm, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;
    double value = -1.;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlange_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlange_Tile_Async( norm, A, &value, sequence, &request );

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
 * @brief Non-blocking equivalent of CHAMELEON_zlange_Tile().
 *
 * @Warning Note that this algorithm includes a RUNTIME_Sequence_wait to cleanup
 * workspaces and thus, will return only when the result is computed. It does
 * not allow to do pipelining.
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
 * @sa CHAMELEON_zlange
 * @sa CHAMELEON_zlange_Tile
 * @sa CHAMELEON_clange_Tile_Async
 * @sa CHAMELEON_dlange_Tile_Async
 * @sa CHAMELEON_slange_Tile_Async
 *
 */
int CHAMELEON_zlange_Tile_Async( cham_normtype_t norm, CHAM_desc_t *A, double *value,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlange_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlange_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlange_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zlange_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zlange_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (norm != ChamMaxNorm) &&
         (norm != ChamOneNorm) &&
         (norm != ChamInfNorm) &&
         (norm != ChamFrobeniusNorm) &&
         (norm != ChamTwoNorm) )
    {
        chameleon_error("CHAMELEON_zlange", "illegal value of norm");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        *value = 0.0;
        return CHAMELEON_SUCCESS;
    }

    if ( norm == ChamTwoNorm ) {
        chameleon_pzgenm2( 1.e-1, A, value, sequence, request );
    }
    else {
        chameleon_pzlange_generic( norm, ChamUpperLower, ChamNonUnit, A, value, sequence, request );
    }

    return CHAMELEON_SUCCESS;
}
