/**
 *
 * @file hip_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_zherk GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "gpuhipblas.h"

int
HIP_zherk( cham_uplo_t uplo, cham_trans_t trans,
           int n, int k,
           const double *alpha,
           const hipDoubleComplex *A, int lda,
           const double *beta,
           hipDoubleComplex *B, int ldb,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasZherk( handle,
                       chameleon_hipblas_const(uplo), chameleon_hipblas_const(trans),
                       n, k,
                       HIPBLAS_VALUE(alpha), A, lda,
                       HIPBLAS_VALUE(beta),  B, ldb );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
