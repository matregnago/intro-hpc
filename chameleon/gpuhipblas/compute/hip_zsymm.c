/**
 *
 * @file hip_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_zsymm GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpuhipblas.h"

int
HIP_zsymm( cham_side_t side, cham_uplo_t uplo,
            int m, int n,
            const hipDoubleComplex *alpha,
            const hipDoubleComplex *A, int lda,
            const hipDoubleComplex *B, int ldb,
            const hipDoubleComplex *beta,
            hipDoubleComplex *C, int ldc,
            hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasZsymm( handle,
                      chameleon_hipblas_const(side), chameleon_hipblas_const(uplo),
                      m, n,
                      HIPBLAS_VALUE(alpha), A, lda,
                                           B, ldb,
                      HIPBLAS_VALUE(beta),  C, ldc );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
