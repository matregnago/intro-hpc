/**
 *
 * @file gpuhipblas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpuhipblas_z_h_
#define _gpuhipblas_z_h_

/**
 *  Declarations of hip kernels - alphabetical order
 */
int HIP_zgemm(  cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *B, int ldb, const hipDoubleComplex *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zhemm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *B, int ldb, const hipDoubleComplex *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zher2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *B, int ldb, const double *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zherk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const hipDoubleComplex *A, int lda, const double *beta, hipDoubleComplex *B, int ldb, hipblasHandle_t handle );
int HIP_zsymm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *B, int ldb, const hipDoubleComplex *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *B, int ldb, const hipDoubleComplex *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zsyrk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, const hipDoubleComplex *beta, hipDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_ztrmm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, hipDoubleComplex *B, int ldb, hipblasHandle_t handle );
int HIP_ztrsm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const hipDoubleComplex *alpha, const hipDoubleComplex *A, int lda, hipDoubleComplex *B, int ldb, hipblasHandle_t handle );

/* hipsolver */
int HIP_zpotrf( cham_uplo_t uplo, int n, hipDoubleComplex *dA, int lda, hipDoubleComplex *dW, int lwork, int *d_info, hipsolverDnHandle_t handle );

#endif /* _gpuhipblas_z_h_ */
