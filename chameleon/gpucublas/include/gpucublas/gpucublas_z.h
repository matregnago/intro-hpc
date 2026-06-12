/**
 *
 * @file gpucublas_z.h
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
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpucublas_z_h_
#define _gpucublas_z_h_

/**
 *  Declarations of cuda kernels - alphabetical order
 */
#if defined(GPUCUBLAS_HAVE_CUDA_HALF)
#if defined(PRECISION_d) || defined(PRECISION_s)
int CUDA_dlag2h( int m, int n, const double *A, int lda, CHAMELEON_Real16_t *B, int ldb, cublasHandle_t handle );
int CUDA_hlag2d( int m, int n, const CHAMELEON_Real16_t *A, int lda, double *B, int ldb, cublasHandle_t handle );
#endif
#endif
int CUDA_zgeadd( cham_trans_t trans, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *beta, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
int CUDA_zgemerge( cham_side_t side, cham_diag_t diag, int M, int N, const cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cublasHandle_t handle );
int CUDA_zgemm(  cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
#if defined(PRECISION_z) || defined(PRECISION_c)
int CUDA_zhemm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
int CUDA_zher2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, const double *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
int CUDA_zherk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const cuDoubleComplex *A, int lda, const double *beta, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
#endif
int CUDA_zherfb( cham_uplo_t uplo, int n, int k, int ib, int nb, const cuDoubleComplex *A, int lda, const cuDoubleComplex *T, int ldt, cuDoubleComplex *C, int ldc, cuDoubleComplex *WORK, int ldwork, cublasHandle_t handle );
int CUDA_zlarfb( cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *C, int LDC, cuDoubleComplex *WORK, int LDWORK, cublasHandle_t handle );
int CUDA_zlacpy( cham_uplo_t uplo, int M, int N, const cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cublasHandle_t handle );
int CUDA_zlaset( cham_uplo_t uplo, int n1, int n2, const cuDoubleComplex *alpha, const cuDoubleComplex *beta, cuDoubleComplex *tileA, int ldtilea, cublasHandle_t handle );
int CUDA_zlatro( cham_uplo_t uplo, cham_trans_t trans, int m, int n, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
int CUDA_zparfb( cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, cublasHandle_t handle );
int CUDA_zpotrf( cham_uplo_t uplo, int n, cuDoubleComplex *A, int lda, cuDoubleComplex *WORK, int lwork, int *d_info, cusolverDnHandle_t handle);
int CUDA_zsymm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
int CUDA_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
int CUDA_zsyrk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, cublasHandle_t handle );
int CUDA_ztpmqrt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cuDoubleComplex *WORK, int lwork, cublasHandle_t handle );
int CUDA_ztpmlqt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cuDoubleComplex *WORK, int lwork, cublasHandle_t handle );
int CUDA_ztrmm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
int CUDA_ztrsm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
int CUDA_ztsmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, cublasHandle_t handle );
int CUDA_zttmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, cublasHandle_t handle );
int CUDA_ztsmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, cublasHandle_t handle );
int CUDA_zttmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, cublasHandle_t handle );
int CUDA_zunmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, cublasHandle_t handle );
int CUDA_zunmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, cublasHandle_t handle );

#endif /* _gpucublas_z_h_ */
