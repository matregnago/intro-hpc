/**
 *
 * @file cuda_zpotrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zpotrf GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

int CUDA_zpotrf( cham_uplo_t uplo, int n,
                 cuDoubleComplex *dA, int lda,
                 cuDoubleComplex *dW, int lwork,
                 int *d_info,
                 cusolverDnHandle_t handle )
{
    cusolverStatus_t rc;

#if !defined(NDEBUG)
    {
        int query;
        rc = cusolverDnZpotrf_bufferSize( handle,
                                          chameleon_cublas_const(uplo),
                                          n, dA, lda, &query );
        assert( rc == CUSOLVER_STATUS_SUCCESS );
        assert( query <= lwork );
    }
#endif

    rc = cusolverDnZpotrf( handle,
                           chameleon_cublas_const(uplo),
                           n, dA, lda, dW, lwork,
                           d_info );

    assert( rc == CUSOLVER_STATUS_SUCCESS );

    return (rc == CUSOLVER_STATUS_SUCCESS) ? CHAMELEON_SUCCESS : CHAMELEON_ERR_UNEXPECTED;
}
