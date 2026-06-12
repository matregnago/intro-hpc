/**
 *
 * @file hip_hgemm.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_hgemm GPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include "gpuhipblas.h"

int
HIP_hgemm( cham_trans_t transa, cham_trans_t transb,
           int m, int n, int k,
           const hipblasHalf *alpha,
           const hipblasHalf *A, int lda,
           const hipblasHalf *B, int ldb,
           const hipblasHalf *beta,
           hipblasHalf *C, int ldc,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasHgemm( handle,
                       chameleon_hipblas_const(transa), chameleon_hipblas_const(transb),
                       m, n, k,
                       HIPBLAS_VALUE(alpha), A, lda,
                                             B, ldb,
                       HIPBLAS_VALUE(beta),  C, ldc );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
