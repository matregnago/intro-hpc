/**
 *
 * @file cuda_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zhemm GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "gpucublas.h"

int
CUDA_zhemm( cham_side_t side, cham_uplo_t uplo,
            int m, int n,
            const cuDoubleComplex *alpha,
            const cuDoubleComplex *A, int lda,
            const cuDoubleComplex *B, int ldb,
            const cuDoubleComplex *beta,
            cuDoubleComplex *C, int ldc,
            cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasZhemm( handle,
                      chameleon_cublas_const(side), chameleon_cublas_const(uplo),
                      m, n,
                      CUBLAS_VALUE(alpha), A, lda,
                                           B, ldb,
                      CUBLAS_VALUE(beta),  C, ldc );

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
