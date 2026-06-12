/**
 *
 * @file cuda_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgemm GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

int
CUDA_zgemm( cham_trans_t transa, cham_trans_t transb,
            int m, int n, int k,
            const cuDoubleComplex *alpha,
            const cuDoubleComplex *A, int lda,
            const cuDoubleComplex *B, int ldb,
            const cuDoubleComplex *beta,
            cuDoubleComplex *C, int ldc,
            cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasZgemm( handle,
                      chameleon_cublas_const(transa), chameleon_cublas_const(transb),
                      m, n, k,
                      CUBLAS_VALUE(alpha), A, lda,
                                           B, ldb,
                      CUBLAS_VALUE(beta),  C, ldc);

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
