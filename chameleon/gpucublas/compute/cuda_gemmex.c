/**
 *
 * @file cuda_gemmex.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_gemmex GPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-07-15
 *
 */
#include "gpucublas.h"

#if !defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
#error "This file should not be compiled"
#endif

int
CUDA_gemmex( cham_trans_t transa, cham_trans_t transb,
             int m, int n, int k,
             const void *alpha,
             const void *A, int lda, cham_flttype_t Atype,
             const void *B, int ldb, cham_flttype_t Btype,
             const void *beta,
             void *C, int ldc, cham_flttype_t Ctype,
             cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasGemmEx( handle,
                       chameleon_cublas_const(transa), chameleon_cublas_const(transb),
                       m, n, k,
                       CUBLAS_VALUE(alpha), A, lda, chameleon_cublas_dtype( Atype ),
                                            B, ldb, chameleon_cublas_dtype( Btype ),
                       CUBLAS_VALUE(beta),  C, ldc, chameleon_cublas_dtype( Ctype ),
                       chameleon_cublas_ctype( Ctype ),
                       CUBLAS_GEMM_DEFAULT );

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
