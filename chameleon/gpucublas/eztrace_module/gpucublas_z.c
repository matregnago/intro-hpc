/**
 *
 * @file gpucublas_z.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Eztrace module file generated from gpucublas_z.h
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> z c d s
 *
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>
#include "gpucublas.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _gpucublas_z_initialized = 0;

#define CURRENT_MODULE gpucublas_z
DECLARE_CURRENT_MODULE;

/* pointers to actual gpucublas_z functions */

#if defined(GPUCUBLAS_HAVE_CUDA_HALF)
#if defined(PRECISION_d) || defined(PRECISION_s)
int (*libCUDA_dlag2h)(int m, int n, const double *A, int lda, CHAMELEON_Real16_t * B, int ldb, cublasHandle_t handle) = NULL;
int (*libCUDA_hlag2d)(int m, int n, const CHAMELEON_Real16_t * A, int lda, double *B, int ldb, cublasHandle_t handle) = NULL;
#endif
#endif
int (*libCUDA_zgeadd)(cham_trans_t trans, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * beta, cuDoubleComplex * B, int ldb, cublasHandle_t handle) = NULL;
int (*libCUDA_zgemerge)(cham_side_t side, cham_diag_t diag, int M, int N, const cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cublasHandle_t handle) = NULL;
int (*libCUDA_zgemm)(cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
int (*libCUDA_zhemm)(cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
int (*libCUDA_zher2k)(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const double *beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
int (*libCUDA_zherk)(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const cuDoubleComplex * A, int lda, const double *beta, cuDoubleComplex * B, int ldb, cublasHandle_t handle) = NULL;
#endif
int (*libCUDA_zherfb)(cham_uplo_t uplo, int n, int k, int ib, int nb, const cuDoubleComplex * A, int lda, const cuDoubleComplex * T, int ldt, cuDoubleComplex * C, int ldc, cuDoubleComplex * WORK, int ldwork, cublasHandle_t handle) = NULL;
int (*libCUDA_zlarfb)(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zparfb)(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zpotrf)(cham_uplo_t uplo, int n, cuDoubleComplex * A, int lda, cuDoubleComplex * WORK, int LWORK, int *d_info, cusolverDnHandle_t handle) = NULL;
int (*libCUDA_zsymm)(cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
int (*libCUDA_zsyr2k)(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
int (*libCUDA_zsyrk)(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) = NULL;
int (*libCUDA_ztpmqrt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cuDoubleComplex * WORK, int lwork, cublasHandle_t handle) = NULL;
int (*libCUDA_ztpmlqt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cuDoubleComplex * WORK, int lwork, cublasHandle_t handle) = NULL;
int (*libCUDA_ztrmm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) = NULL;
int (*libCUDA_ztrsm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) = NULL;
int (*libCUDA_ztsmlq)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zttmlq)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_ztsmqr)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zttmqr)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zunmlqt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex * A, int LDA, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) = NULL;
int (*libCUDA_zunmqrt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex * A, int LDA, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) = NULL;



/* Wrapper functions */

#if defined(GPUCUBLAS_HAVE_CUDA_HALF)
#if defined(PRECISION_d) || defined(PRECISION_s)

int CUDA_dlag2h(int m, int n, const double *A, int lda, CHAMELEON_Real16_t * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_dlag2h(m, n, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_hlag2d(int m, int n, const CHAMELEON_Real16_t * A, int lda, double *B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_hlag2d(m, n, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}
#endif
#endif

int CUDA_zgeadd(cham_trans_t trans, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * beta, cuDoubleComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zgeadd(trans, m, n, alpha, A, lda, beta, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zgemerge(cham_side_t side, cham_diag_t diag, int M, int N, const cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zgemerge(side, diag, M, N, A, LDA, B, LDB, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zgemm(cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}
#if defined(PRECISION_z) || defined(PRECISION_c)

int CUDA_zhemm(cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zher2k(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const double *beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zher2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zherk(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const cuDoubleComplex * A, int lda, const double *beta, cuDoubleComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zherk(uplo, trans, n, k, alpha, A, lda, beta, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}
#endif

int CUDA_zherfb(cham_uplo_t uplo, int n, int k, int ib, int nb, const cuDoubleComplex * A, int lda, const cuDoubleComplex * T, int ldt, cuDoubleComplex * C, int ldc, cuDoubleComplex * WORK, int ldwork, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zherfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zlarfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zlarfb(side, trans, direct, storev, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zparfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zparfb(side, trans, direct, storev, M1, N1, M2, N2, K, L, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zpotrf(cham_uplo_t uplo, int n, cuDoubleComplex *A, int lda, cuDoubleComplex *WORK, int lwork, int *d_info, cusolverDnHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zpotrf(uplo, n, A, lda, WORK, lwork, d_info, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zsymm(cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zsyr2k(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * B, int ldb, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zsyrk(cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, const cuDoubleComplex * beta, cuDoubleComplex * C, int ldc, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztpmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cuDoubleComplex * WORK, int lwork, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztpmqrt(side, trans, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, lwork, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztpmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * A, int LDA, cuDoubleComplex * B, int LDB, cuDoubleComplex * WORK, int lwork, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztpmlqt(side, trans, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, lwork, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const cuDoubleComplex * alpha, const cuDoubleComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztsmlq(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztsmlq(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zttmlq(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zttmlq(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_ztsmqr(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_ztsmqr(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zttmqr(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex * A1, int LDA1, cuDoubleComplex * A2, int LDA2, const cuDoubleComplex * V, int LDV, const cuDoubleComplex * T, int LDT, cuDoubleComplex * WORK, int LWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zttmqr(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zunmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex * A, int LDA, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zunmlqt(side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zunmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex * A, int LDA, const cuDoubleComplex * T, int LDT, cuDoubleComplex * C, int LDC, cuDoubleComplex * WORK, int LDWORK, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zunmqrt(side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK, handle);
    FUNCTION_EXIT;
    return ret;
}


PPTRACE_START_INTERCEPT_FUNCTIONS(gpucublas_z)
#if defined(GPUCUBLAS_HAVE_CUDA_HALF)
#if defined(PRECISION_d) || defined(PRECISION_s)
INTERCEPT3("CUDA_dlag2h", libCUDA_dlag2h)
INTERCEPT3("CUDA_hlag2d", libCUDA_hlag2d)
#endif
#endif
INTERCEPT3("CUDA_zgeadd", libCUDA_zgeadd)
INTERCEPT3("CUDA_zgemerge", libCUDA_zgemerge)
INTERCEPT3("CUDA_zgemm", libCUDA_zgemm)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("CUDA_zhemm", libCUDA_zhemm)
INTERCEPT3("CUDA_zher2k", libCUDA_zher2k)
INTERCEPT3("CUDA_zherk", libCUDA_zherk)
#endif
INTERCEPT3("CUDA_zherfb", libCUDA_zherfb)
INTERCEPT3("CUDA_zlarfb", libCUDA_zlarfb)
INTERCEPT3("CUDA_zparfb", libCUDA_zparfb)
INTERCEPT3("CUDA_zsymm", libCUDA_zsymm)
INTERCEPT3("CUDA_zsyr2k", libCUDA_zsyr2k)
INTERCEPT3("CUDA_zsyrk", libCUDA_zsyrk)
INTERCEPT3("CUDA_ztpmqrt", libCUDA_ztpmqrt)
INTERCEPT3("CUDA_ztpmlqt", libCUDA_ztpmlqt)
INTERCEPT3("CUDA_ztrmm", libCUDA_ztrmm)
INTERCEPT3("CUDA_ztrsm", libCUDA_ztrsm)
INTERCEPT3("CUDA_ztsmlq", libCUDA_ztsmlq)
INTERCEPT3("CUDA_zttmlq", libCUDA_zttmlq)
INTERCEPT3("CUDA_ztsmqr", libCUDA_ztsmqr)
INTERCEPT3("CUDA_zttmqr", libCUDA_zttmqr)
INTERCEPT3("CUDA_zunmlqt", libCUDA_zunmlqt)
INTERCEPT3("CUDA_zunmqrt", libCUDA_zunmqrt)

PPTRACE_END_INTERCEPT_FUNCTIONS(gpucublas_z)
