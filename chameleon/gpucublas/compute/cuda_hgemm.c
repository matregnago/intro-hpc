/**
 *
 * @file cuda_hgemm.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_hgemm GPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-07-15
 *
 */
#include "gpucublas.h"

#if !defined(GPUCUBLAS_HAVE_CUBLASHGEMM)
#error "This file should not be compiled"
#endif

extern "C" int
CUDA_hgemm( cham_trans_t transa, cham_trans_t transb,
            int m, int n, int k,
            const CHAMELEON_Real16_t *alpha,
            const CHAMELEON_Real16_t *A, int lda,
            const CHAMELEON_Real16_t *B, int ldb,
            const CHAMELEON_Real16_t *beta,
            CHAMELEON_Real16_t *C, int ldc,
            cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasHgemm( handle,
                      (cublasOperation_t)chameleon_cublas_const(transa),
                      (cublasOperation_t)chameleon_cublas_const(transb),
                      m, n, k,
                      CUBLAS_VALUE(alpha), A, lda,
                                           B, ldb,
                      CUBLAS_VALUE(beta),  C, ldc);

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
