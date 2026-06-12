/**
 *
 * @file hip_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_ztrmm GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @author Matthieu Kuhn
 * @author Xavier Lacoste
 * @date 2024-04-24
 * @precisions normal z -> c d s
 *
 */
#include "gpuhipblas.h"

int
HIP_ztrmm( cham_side_t side, cham_uplo_t uplo,
           cham_trans_t transa, cham_diag_t diag,
           int m, int n,
           const hipDoubleComplex *alpha,
           const hipDoubleComplex *A, int lda,
           hipDoubleComplex *B, int ldb,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

#if (HIP_VERSION_MAJOR < 6)
    rc = hipblasZtrmm( handle,
                       chameleon_hipblas_const(side), chameleon_hipblas_const(uplo),
                       chameleon_hipblas_const(transa), chameleon_hipblas_const(diag),
                       m, n,
                       HIPBLAS_VALUE(alpha), A, lda,
                       B, ldb );
#else
    rc = hipblasZtrmm( handle,
                       chameleon_hipblas_const(side), chameleon_hipblas_const(uplo),
                       chameleon_hipblas_const(transa), chameleon_hipblas_const(diag),
                       m, n,
                       HIPBLAS_VALUE(alpha), A, lda,
                       B, ldb,
                       B, ldb );
#endif

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
