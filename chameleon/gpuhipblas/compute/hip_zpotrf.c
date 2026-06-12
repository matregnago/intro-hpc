/**
 *
 * @file hip_zpotrf.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_zpotrf GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2026-02-16
 * @precisions normal z -> c d s
 *
 */
#include "gpuhipblas.h"

int
HIP_zpotrf( cham_uplo_t uplo, int n,
            hipDoubleComplex *dA, int lda,
            hipDoubleComplex *dW, int lwork,
            int *d_info,
            hipsolverDnHandle_t handle )
{
    hipsolverStatus_t rc;

#if !defined(NDEBUG)
    {
        int query;
        rc = hipsolverDnZpotrf_bufferSize( handle,
                                           chameleon_hipblas_const(uplo),
                                           n, dA, lda, &query );
        assert( rc == HIPSOLVER_STATUS_SUCCESS );
        assert( query <= lwork );
    }
#endif

    rc = hipsolverDnZpotrf( handle,
                            chameleon_hipblas_const(uplo),
                            n, dA, lda, dW, lwork,
                            d_info );

    assert( rc == HIPSOLVER_STATUS_SUCCESS );

    return (rc == HIPSOLVER_STATUS_SUCCESS) ? CHAMELEON_SUCCESS : CHAMELEON_ERR_UNEXPECTED;
}
