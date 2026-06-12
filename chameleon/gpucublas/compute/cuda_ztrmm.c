/**
 *
 * @file cuda_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztrmm GPU kernel
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
CUDA_ztrmm( cham_side_t side, cham_uplo_t uplo,
            cham_trans_t transa, cham_diag_t diag,
            int m, int n,
            const cuDoubleComplex *alpha,
            const cuDoubleComplex *A, int lda,
            cuDoubleComplex *B, int ldb,
            cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasZtrmm( handle,
                      chameleon_cublas_const(side), chameleon_cublas_const(uplo),
                      chameleon_cublas_const(transa), chameleon_cublas_const(diag),
                      m, n,
                      CUBLAS_VALUE(alpha), A, lda,
                      B, ldb,
                      B, ldb );

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}

