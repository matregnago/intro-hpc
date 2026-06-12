/**
 *
 * @file tcore_z.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Eztrace module file generated from tcore_z.h
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @author Florent Pruvost
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
#include "coreblas.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _coreblas_tcore_z_initialized = 0;

#define CURRENT_MODULE coreblas_tcore_z
DECLARE_CURRENT_MODULE;

/* pointers to actual coreblas_tcore_z functions */

#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libTCORE_dlag2z)(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
#endif
void (*libTCORE_dzasum)(cham_store_t storev, cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, double *work) = NULL;
int (*libTCORE_zaxpy)(int M, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, int incA, CHAM_tile_t * B, int incB) = NULL;
int (*libTCORE_zgeadd)(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * B) = NULL;
int (*libTCORE_zgelqt)(int M, int N, int IB, CHAM_tile_t * A, CHAM_tile_t * T, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
void (*libTCORE_zgemv)(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * x, int incx, CHAMELEON_Complex64_t beta, CHAM_tile_t * y, int incy) = NULL;
void (*libTCORE_zgemm)(cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) = NULL;
int (*libTCORE_zgeqrt)(int M, int N, int IB, CHAM_tile_t * A, CHAM_tile_t * T, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libTCORE_zgessm)(int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t * L, CHAM_tile_t * A) = NULL;
int (*libTCORE_zgessq)(cham_store_t storev, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) = NULL;
int (*libTCORE_zgetrf)(int M, int N, CHAM_tile_t * A, int *IPIV, int *INFO) = NULL;
int (*libTCORE_zgetrf_incpiv)(int M, int N, int IB, CHAM_tile_t * A, int *IPIV, int *INFO) = NULL;
int (*libTCORE_zgetrf_nopiv)(int M, int N, int IB, CHAM_tile_t * A, int *INFO) = NULL;
void (*libTCORE_zhe2ge)(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libTCORE_zhemm)(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) = NULL;
void (*libTCORE_zherk)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, const CHAM_tile_t * A, double beta, CHAM_tile_t * C) = NULL;
void (*libTCORE_zher2k)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, double beta, CHAM_tile_t * C) = NULL;
#endif
int (*libTCORE_zherfb)(cham_uplo_t uplo, int N, int K, int IB, int NB, const CHAM_tile_t * A, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int ldwork) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
int (*libTCORE_zhessq)(cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) = NULL;
#endif
void (*libTCORE_zlacpy)(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
void (*libTCORE_zlacpyx)(cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t * A, int LDA, int displB, CHAM_tile_t * B, int LDB) = NULL;
void (*libTCORE_zlange)(cham_normtype_t norm, int M, int N, const CHAM_tile_t * A, double *work, double *normA) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libTCORE_zlanhe)(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t * A, double *work, double *normA) = NULL;
#endif
void (*libTCORE_zlansy)(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t * A, double *work, double *normA) = NULL;
void (*libTCORE_zlantr)(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, double *work, double *normA) = NULL;
int (*libTCORE_zlascal)(cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t * A) = NULL;
void (*libTCORE_zlaset)(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_tile_t * A) = NULL;
void (*libTCORE_zlaset2)(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAM_tile_t * A) = NULL;
int (*libTCORE_zlaswp_get)(cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t * A, CHAM_tile_t * B, const int *perm) = NULL;
int (*libTCORE_zlaswp_set)(cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t * A, CHAM_tile_t * B, const int *invp) = NULL;
int (*libTCORE_zlatro)(cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
void (*libTCORE_zlauum)(cham_uplo_t uplo, int N, CHAM_tile_t * A) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libTCORE_zplghe)(double bump, int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
#endif
void (*libTCORE_zplgsy)(CHAMELEON_Complex64_t bump, int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
void (*libTCORE_zplrnt)(int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
void (*libTCORE_zpotrf)(cham_uplo_t uplo, int n, CHAM_tile_t * A, int *INFO) = NULL;
int (*libTCORE_zssssm)(int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * L1, const CHAM_tile_t * L2, const int *IPIV) = NULL;
void (*libTCORE_zsymm)(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) = NULL;
void (*libTCORE_zsyrk)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) = NULL;
void (*libTCORE_zsyr2k)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) = NULL;
int (*libTCORE_zsyssq)(cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) = NULL;
int (*libTCORE_zsytf2_nopiv)(cham_uplo_t uplo, int n, CHAM_tile_t * A) = NULL;
int (*libTCORE_ztplqt)(int M, int N, int L, int IB, CHAM_tile_t * A, CHAM_tile_t * B, CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libTCORE_ztpmlqt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * A, CHAM_tile_t * B, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libTCORE_ztpmqrt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * A, CHAM_tile_t * B, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libTCORE_ztpqrt)(int M, int N, int L, int IB, CHAM_tile_t * A, CHAM_tile_t * B, CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libTCORE_ztradd)(cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * B) = NULL;
void (*libTCORE_ztrasm)(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, double *work) = NULL;
void (*libTCORE_ztrmm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
void (*libTCORE_ztrsm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
int (*libTCORE_ztrssq)(cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) = NULL;
void (*libTCORE_ztrtri)(cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t * A, int *info) = NULL;
int (*libTCORE_ztsmlq_hetra1)(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK, int ldwork) = NULL;
int (*libTCORE_ztsmqr_hetra1)(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK, int ldwork) = NULL;
int (*libTCORE_ztstrf)(int M, int N, int IB, int NB, CHAM_tile_t * U, CHAM_tile_t * A, CHAM_tile_t * L, int *IPIV, CHAMELEON_Complex64_t * WORK, int LDWORK, int *INFO) = NULL;
int (*libTCORE_zunmlq)(cham_side_t side, cham_trans_t trans, int M, int N, int IB, int K, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libTCORE_zunmqr)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libTCORE_zgesum)(cham_store_t storev, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sum) = NULL;
int (*libTCORE_zcesca)(int center, int scale, cham_store_t axis, int M, int N, int Mt, int Nt, const CHAM_tile_t * Gi, const CHAM_tile_t * Gj, const CHAM_tile_t * G, const CHAM_tile_t * Di, const CHAM_tile_t * Dj, CHAM_tile_t * A) = NULL;
int (*libTCORE_zgram)(cham_uplo_t uplo, int M, int N, int Mt, int Nt, const CHAM_tile_t * Di, const CHAM_tile_t * Dj, const CHAM_tile_t * D, CHAM_tile_t * A) = NULL;
void (*libTCORE_zprint)(FILE * file, const char *header, cham_uplo_t uplo, int M, int N, int Am, int An, const CHAM_tile_t * A) = NULL;



/* Wrapper functions */

#if defined(PRECISION_z) || defined(PRECISION_c)

void TCORE_dlag2z(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_dlag2z(uplo, M, N, A, B);
    FUNCTION_EXIT;

}
#endif

void TCORE_dzasum(cham_store_t storev, cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, double *work) {
    FUNCTION_ENTRY;
    libTCORE_dzasum(storev, uplo, M, N, A, work);
    FUNCTION_EXIT;

}

int TCORE_zaxpy(int M, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, int incA, CHAM_tile_t * B, int incB) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zaxpy(M, alpha, A, incA, B, incB);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgeadd(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgeadd(trans, M, N, alpha, A, beta, B);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgelqt(int M, int N, int IB, CHAM_tile_t * A, CHAM_tile_t * T, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgelqt(M, N, IB, A, T, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zgemv(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * x, int incx, CHAMELEON_Complex64_t beta, CHAM_tile_t * y, int incy) {
    FUNCTION_ENTRY;
    libTCORE_zgemv(trans, M, N, alpha, A, x, incx, beta, y, incy);
    FUNCTION_EXIT;

}

void TCORE_zgemm(cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zgemm(transA, transB, M, N, K, alpha, A, B, beta, C);
    FUNCTION_EXIT;

}

int TCORE_zgeqrt(int M, int N, int IB, CHAM_tile_t * A, CHAM_tile_t * T, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgeqrt(M, N, IB, A, T, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgessm(int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t * L, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgessm(M, N, K, IB, IPIV, L, A);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgessq(cham_store_t storev, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgessq(storev, M, N, A, sclssq);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgetrf(int M, int N, CHAM_tile_t * A, int *IPIV, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgetrf(M, N, A, IPIV, INFO);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgetrf_incpiv(int M, int N, int IB, CHAM_tile_t * A, int *IPIV, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgetrf_incpiv(M, N, IB, A, IPIV, INFO);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgetrf_nopiv(int M, int N, int IB, CHAM_tile_t * A, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgetrf_nopiv(M, N, IB, A, INFO);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zhe2ge(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_zhe2ge(uplo, M, N, A, B);
    FUNCTION_EXIT;

}
#if defined(PRECISION_z) || defined(PRECISION_c)

void TCORE_zhemm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zhemm(side, uplo, M, N, alpha, A, B, beta, C);
    FUNCTION_EXIT;

}

void TCORE_zherk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, const CHAM_tile_t * A, double beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zherk(uplo, trans, N, K, alpha, A, beta, C);
    FUNCTION_EXIT;

}

void TCORE_zher2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, double beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zher2k(uplo, trans, N, K, alpha, A, B, beta, C);
    FUNCTION_EXIT;

}
#endif

int TCORE_zherfb(cham_uplo_t uplo, int N, int K, int IB, int NB, const CHAM_tile_t * A, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int ldwork) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zherfb(uplo, N, K, IB, NB, A, T, C, WORK, ldwork);
    FUNCTION_EXIT;
    return ret;
}
#if defined(PRECISION_z) || defined(PRECISION_c)

int TCORE_zhessq(cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zhessq(storev, uplo, N, A, sclssq);
    FUNCTION_EXIT;
    return ret;
}
#endif

void TCORE_zlacpy(cham_uplo_t uplo, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_zlacpy(uplo, M, N, A, B);
    FUNCTION_EXIT;

}

void TCORE_zlacpyx(cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t * A, int LDA, int displB, CHAM_tile_t * B, int LDB) {
    FUNCTION_ENTRY;
    libTCORE_zlacpyx(uplo, M, N, displA, A, LDA, displB, B, LDB);
    FUNCTION_EXIT;

}

void TCORE_zlange(cham_normtype_t norm, int M, int N, const CHAM_tile_t * A, double *work, double *normA) {
    FUNCTION_ENTRY;
    libTCORE_zlange(norm, M, N, A, work, normA);
    FUNCTION_EXIT;

}
#if defined(PRECISION_z) || defined(PRECISION_c)

void TCORE_zlanhe(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t * A, double *work, double *normA) {
    FUNCTION_ENTRY;
    libTCORE_zlanhe(norm, uplo, N, A, work, normA);
    FUNCTION_EXIT;

}
#endif

void TCORE_zlansy(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t * A, double *work, double *normA) {
    FUNCTION_ENTRY;
    libTCORE_zlansy(norm, uplo, N, A, work, normA);
    FUNCTION_EXIT;

}

void TCORE_zlantr(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, double *work, double *normA) {
    FUNCTION_ENTRY;
    libTCORE_zlantr(norm, uplo, diag, M, N, A, work, normA);
    FUNCTION_EXIT;

}

int TCORE_zlascal(cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zlascal(uplo, m, n, alpha, A);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zlaset(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    libTCORE_zlaset(uplo, n1, n2, alpha, beta, A);
    FUNCTION_EXIT;

}

void TCORE_zlaset2(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    libTCORE_zlaset2(uplo, n1, n2, alpha, A);
    FUNCTION_EXIT;

}

int TCORE_zlaswp_get(cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t * A, CHAM_tile_t * B, const int *perm) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zlaswp_get(side, m0, m, n, k, A, B, perm);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zlaswp_set(cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t * A, CHAM_tile_t * B, const int *invp) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zlaswp_set(side, m0, m, n, k, A, B, invp);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zlatro(cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zlatro(uplo, trans, M, N, A, B);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zlauum(cham_uplo_t uplo, int N, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    libTCORE_zlauum(uplo, N, A);
    FUNCTION_EXIT;

}
#if defined(PRECISION_z) || defined(PRECISION_c)

void TCORE_zplghe(double bump, int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libTCORE_zplghe(bump, m, n, tileA, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}
#endif

void TCORE_zplgsy(CHAMELEON_Complex64_t bump, int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libTCORE_zplgsy(bump, m, n, tileA, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}

void TCORE_zplrnt(int m, int n, CHAM_tile_t * tileA, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libTCORE_zplrnt(m, n, tileA, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}

void TCORE_zpotrf(cham_uplo_t uplo, int n, CHAM_tile_t * A, int *INFO) {
    FUNCTION_ENTRY;
    libTCORE_zpotrf(uplo, n, A, INFO);
    FUNCTION_EXIT;

}

int TCORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * L1, const CHAM_tile_t * L2, const int *IPIV) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zssssm(M1, N1, M2, N2, K, IB, A1, A2, L1, L2, IPIV);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zsymm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zsymm(side, uplo, M, N, alpha, A, B, beta, C);
    FUNCTION_EXIT;

}

void TCORE_zsyrk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zsyrk(uplo, trans, N, K, alpha, A, beta, C);
    FUNCTION_EXIT;

}

void TCORE_zsyr2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, const CHAM_tile_t * B, CHAMELEON_Complex64_t beta, CHAM_tile_t * C) {
    FUNCTION_ENTRY;
    libTCORE_zsyr2k(uplo, trans, N, K, alpha, A, B, beta, C);
    FUNCTION_EXIT;

}

int TCORE_zsyssq(cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zsyssq(storev, uplo, N, A, sclssq);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zsytf2_nopiv(cham_uplo_t uplo, int n, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zsytf2_nopiv(uplo, n, A);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztplqt(int M, int N, int L, int IB, CHAM_tile_t * A, CHAM_tile_t * B, CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztplqt(M, N, L, IB, A, B, T, WORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztpmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * A, CHAM_tile_t * B, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztpmlqt(side, trans, M, N, K, L, IB, V, T, A, B, WORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztpmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * A, CHAM_tile_t * B, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztpmqrt(side, trans, M, N, K, L, IB, V, T, A, B, WORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztpqrt(int M, int N, int L, int IB, CHAM_tile_t * A, CHAM_tile_t * B, CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztpqrt(M, N, L, IB, A, B, T, WORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztradd(cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAMELEON_Complex64_t beta, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztradd(uplo, trans, M, N, alpha, A, beta, B);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_ztrasm(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, double *work) {
    FUNCTION_ENTRY;
    libTCORE_ztrasm(storev, uplo, diag, M, N, A, work);
    FUNCTION_EXIT;

}

void TCORE_ztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_ztrmm(side, uplo, transA, diag, M, N, alpha, A, B);
    FUNCTION_EXIT;

}

void TCORE_ztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_ztrsm(side, uplo, transA, diag, M, N, alpha, A, B);
    FUNCTION_EXIT;

}

int TCORE_ztrssq(cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sclssq) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztrssq(uplo, diag, M, N, A, sclssq);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_ztrtri(cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t * A, int *info) {
    FUNCTION_ENTRY;
    libTCORE_ztrtri(uplo, diag, N, A, info);
    FUNCTION_EXIT;

}

int TCORE_ztsmlq_hetra1(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK, int ldwork) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztsmlq_hetra1(side, trans, m1, n1, m2, n2, k, ib, A1, A2, V, T, WORK, ldwork);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztsmqr_hetra1(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t * A1, CHAM_tile_t * A2, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAMELEON_Complex64_t * WORK, int ldwork) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztsmqr_hetra1(side, trans, m1, n1, m2, n2, k, ib, A1, A2, V, T, WORK, ldwork);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_ztstrf(int M, int N, int IB, int NB, CHAM_tile_t * U, CHAM_tile_t * A, CHAM_tile_t * L, int *IPIV, CHAMELEON_Complex64_t * WORK, int LDWORK, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libTCORE_ztstrf(M, N, IB, NB, U, A, L, IPIV, WORK, LDWORK, INFO);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zunmlq(cham_side_t side, cham_trans_t trans, int M, int N, int IB, int K, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zunmlq(side, trans, M, N, IB, K, V, T, C, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zunmqr(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAM_tile_t * V, const CHAM_tile_t * T, CHAM_tile_t * C, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zunmqr(side, trans, M, N, K, IB, V, T, C, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgesum(cham_store_t storev, int M, int N, const CHAM_tile_t * A, CHAM_tile_t * sum) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgesum(storev, M, N, A, sum);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zcesca(int center, int scale, cham_store_t axis, int M, int N, int Mt, int Nt, const CHAM_tile_t * Gi, const CHAM_tile_t * Gj, const CHAM_tile_t * G, const CHAM_tile_t * Di, const CHAM_tile_t * Dj, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zcesca(center, scale, axis, M, N, Mt, Nt, Gi, Gj, G, Di, Dj, A);
    FUNCTION_EXIT;
    return ret;
}

int TCORE_zgram(cham_uplo_t uplo, int M, int N, int Mt, int Nt, const CHAM_tile_t * Di, const CHAM_tile_t * Dj, const CHAM_tile_t * D, CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    int ret = libTCORE_zgram(uplo, M, N, Mt, Nt, Di, Dj, D, A);
    FUNCTION_EXIT;
    return ret;
}

void TCORE_zprint(FILE * file, const char *header, cham_uplo_t uplo, int M, int N, int Am, int An, const CHAM_tile_t * A) {
    FUNCTION_ENTRY;
    libTCORE_zprint(file, header, uplo, M, N, Am, An, A);
    FUNCTION_EXIT;

}


PPTRACE_START_INTERCEPT_FUNCTIONS(coreblas_tcore_z)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("TCORE_dlag2z", libTCORE_dlag2z)
#endif
INTERCEPT3("TCORE_dzasum", libTCORE_dzasum)
INTERCEPT3("TCORE_zaxpy", libTCORE_zaxpy)
INTERCEPT3("TCORE_zgeadd", libTCORE_zgeadd)
INTERCEPT3("TCORE_zgelqt", libTCORE_zgelqt)
INTERCEPT3("TCORE_zgemv", libTCORE_zgemv)
INTERCEPT3("TCORE_zgemm", libTCORE_zgemm)
INTERCEPT3("TCORE_zgeqrt", libTCORE_zgeqrt)
INTERCEPT3("TCORE_zgessm", libTCORE_zgessm)
INTERCEPT3("TCORE_zgessq", libTCORE_zgessq)
INTERCEPT3("TCORE_zgetrf", libTCORE_zgetrf)
INTERCEPT3("TCORE_zgetrf_incpiv", libTCORE_zgetrf_incpiv)
INTERCEPT3("TCORE_zgetrf_nopiv", libTCORE_zgetrf_nopiv)
INTERCEPT3("TCORE_zhe2ge", libTCORE_zhe2ge)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("TCORE_zhemm", libTCORE_zhemm)
INTERCEPT3("TCORE_zherk", libTCORE_zherk)
INTERCEPT3("TCORE_zher2k", libTCORE_zher2k)
#endif
INTERCEPT3("TCORE_zherfb", libTCORE_zherfb)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("TCORE_zhessq", libTCORE_zhessq)
#endif
INTERCEPT3("TCORE_zlacpy", libTCORE_zlacpy)
INTERCEPT3("TCORE_zlacpyx", libTCORE_zlacpyx)
INTERCEPT3("TCORE_zlange", libTCORE_zlange)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("TCORE_zlanhe", libTCORE_zlanhe)
#endif
INTERCEPT3("TCORE_zlansy", libTCORE_zlansy)
INTERCEPT3("TCORE_zlantr", libTCORE_zlantr)
INTERCEPT3("TCORE_zlascal", libTCORE_zlascal)
INTERCEPT3("TCORE_zlaset", libTCORE_zlaset)
INTERCEPT3("TCORE_zlaset2", libTCORE_zlaset2)
INTERCEPT3("TCORE_zlaswp_get", libTCORE_zlaswp_get)
INTERCEPT3("TCORE_zlaswp_set", libTCORE_zlaswp_set)
INTERCEPT3("TCORE_zlatro", libTCORE_zlatro)
INTERCEPT3("TCORE_zlauum", libTCORE_zlauum)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("TCORE_zplghe", libTCORE_zplghe)
#endif
INTERCEPT3("TCORE_zplgsy", libTCORE_zplgsy)
INTERCEPT3("TCORE_zplrnt", libTCORE_zplrnt)
INTERCEPT3("TCORE_zpotrf", libTCORE_zpotrf)
INTERCEPT3("TCORE_zssssm", libTCORE_zssssm)
INTERCEPT3("TCORE_zsymm", libTCORE_zsymm)
INTERCEPT3("TCORE_zsyrk", libTCORE_zsyrk)
INTERCEPT3("TCORE_zsyr2k", libTCORE_zsyr2k)
INTERCEPT3("TCORE_zsyssq", libTCORE_zsyssq)
INTERCEPT3("TCORE_zsytf2_nopiv", libTCORE_zsytf2_nopiv)
INTERCEPT3("TCORE_ztplqt", libTCORE_ztplqt)
INTERCEPT3("TCORE_ztpmlqt", libTCORE_ztpmlqt)
INTERCEPT3("TCORE_ztpmqrt", libTCORE_ztpmqrt)
INTERCEPT3("TCORE_ztpqrt", libTCORE_ztpqrt)
INTERCEPT3("TCORE_ztradd", libTCORE_ztradd)
INTERCEPT3("TCORE_ztrasm", libTCORE_ztrasm)
INTERCEPT3("TCORE_ztrmm", libTCORE_ztrmm)
INTERCEPT3("TCORE_ztrsm", libTCORE_ztrsm)
INTERCEPT3("TCORE_ztrssq", libTCORE_ztrssq)
INTERCEPT3("TCORE_ztrtri", libTCORE_ztrtri)
INTERCEPT3("TCORE_ztsmlq_hetra1", libTCORE_ztsmlq_hetra1)
INTERCEPT3("TCORE_ztsmqr_hetra1", libTCORE_ztsmqr_hetra1)
INTERCEPT3("TCORE_ztstrf", libTCORE_ztstrf)
INTERCEPT3("TCORE_zunmlq", libTCORE_zunmlq)
INTERCEPT3("TCORE_zunmqr", libTCORE_zunmqr)
INTERCEPT3("TCORE_zgesum", libTCORE_zgesum)
INTERCEPT3("TCORE_zcesca", libTCORE_zcesca)
INTERCEPT3("TCORE_zgram", libTCORE_zgram)
INTERCEPT3("TCORE_zprint", libTCORE_zprint)

PPTRACE_END_INTERCEPT_FUNCTIONS(coreblas_tcore_z)
