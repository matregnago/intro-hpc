/**
 *
 * @file core_z.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Eztrace module file generated from core_z.h
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
static volatile int _coreblas_core_z_initialized = 0;

#define CURRENT_MODULE coreblas_core_z
DECLARE_CURRENT_MODULE;

/* pointers to actual coreblas_core_z functions */

void (*libCORE_dlag2z)(cham_uplo_t uplo, int M, int N, const double *A, int lda, CHAMELEON_Complex64_t * B, int ldb) = NULL;
void (*libCORE_dzasum)(cham_store_t storev, cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int lda, double *work) = NULL;
int (*libCORE_zaxpy)(int M, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int incA, CHAMELEON_Complex64_t * B, int incB) = NULL;
int (*libCORE_zgeadd)(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * B, int LDB) = NULL;
int (*libCORE_zlascal)(cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t * A, int lda) = NULL;
int (*libCORE_zgelqt)(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
void (*libCORE_zgemm)(cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
void (*libCORE_zgemv)(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * x, int incX, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * y, int incY) = NULL;
int (*libCORE_zgeqrt)(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_zgesplit)(cham_side_t side, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
int (*libCORE_zgessm)(int M, int N, int K, int IB, const int *IPIV, const CHAMELEON_Complex64_t * L, int LDL, CHAMELEON_Complex64_t * A, int LDA) = NULL;
int (*libCORE_zgessq)(cham_store_t storev, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) = NULL;
int (*libCORE_zgetf2_nopiv)(int M, int N, CHAMELEON_Complex64_t * A, int LDA) = NULL;
int (*libCORE_zgetrf)(int M, int N, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *INFO) = NULL;
int (*libCORE_zgetrf_incpiv)(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *INFO) = NULL;
int (*libCORE_zgetrf_nopiv)(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, int *INFO) = NULL;
int (*libCORE_zgetrf_panel_diag)(int m, int n, int h, int m0, int ib, CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * U, int ldu, int *IPIV, CHAM_pivot_t * nextpiv, const CHAM_pivot_t * prevpiv) = NULL;
int (*libCORE_zgetrf_panel_offdiag)(int m, int n, int h, int m0, int ib, CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * U, int ldu, CHAM_pivot_t * nextpiv, const CHAM_pivot_t * prevpiv) = NULL;
int (*libCORE_zgetrf_reclap)(int M, int N, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *info) = NULL;
int (*libCORE_zgetrf_rectil)(const CHAM_desc_t A, int *IPIV, int *info) = NULL;
void (*libCORE_zgetrip)(int m, int n, CHAMELEON_Complex64_t * A, CHAMELEON_Complex64_t * work) = NULL;
void (*libCORE_zhe2ge)(cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
void (*libCORE_zhegst)(int itype, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, int *INFO) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libCORE_zhemm)(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
void (*libCORE_zherk)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, const CHAMELEON_Complex64_t * A, int LDA, double beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
void (*libCORE_zher2k)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, double beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
int (*libCORE_zhessq)(cham_store_t storev, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) = NULL;
#endif
int (*libCORE_zherfb)(cham_uplo_t uplo, int N, int K, int IB, int NB, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
void (*libCORE_zlacpy)(cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
void (*libCORE_zlange)(cham_normtype_t norm, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libCORE_zlanhe)(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) = NULL;
#endif
void (*libCORE_zlansy)(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) = NULL;
void (*libCORE_zlantr)(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) = NULL;
int (*libCORE_zlarfb_gemm)(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zlarfx2)(cham_side_t side, int N, CHAMELEON_Complex64_t V, CHAMELEON_Complex64_t TAU, CHAMELEON_Complex64_t * C1, int LDC1, CHAMELEON_Complex64_t * C2, int LDC2) = NULL;
int (*libCORE_zlarfx2c)(cham_uplo_t uplo, CHAMELEON_Complex64_t V, CHAMELEON_Complex64_t TAU, CHAMELEON_Complex64_t * C1, CHAMELEON_Complex64_t * C2, CHAMELEON_Complex64_t * C3) = NULL;
int (*libCORE_zlarfx2ce)(cham_uplo_t uplo, CHAMELEON_Complex64_t * V, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * C1, CHAMELEON_Complex64_t * C2, CHAMELEON_Complex64_t * C3) = NULL;
void (*libCORE_zlarfy)(int N, CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * V, const CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
void (*libCORE_zlaset)(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * tileA, int ldtilea) = NULL;
void (*libCORE_zlaset2)(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t * tileA, int ldtilea) = NULL;
void (*libCORE_zlaswp)(int N, CHAMELEON_Complex64_t * A, int LDA, int I1, int I2, const int *IPIV, int INC) = NULL;
int (*libCORE_zlaswp_get)(cham_side_t side, int m0, int m, int n, int k, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb, const int *perm) = NULL;
int (*libCORE_zlaswp_set)(cham_side_t side, int m0, int m, int n, int k, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb, const int *invp) = NULL;
int (*libCORE_zlaswp_ontile)(CHAM_desc_t descA, int i1, int i2, const int *ipiv, int inc) = NULL;
int (*libCORE_zlatro)(cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
int (*libCORE_zlatm1)(int MODE, double COND, int IRSIGN, cham_dist_t DIST, unsigned long long int seed, CHAMELEON_Complex64_t * D, int N) = NULL;
void (*libCORE_zlauum)(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA) = NULL;
int (*libCORE_zpamm)(int op, cham_side_t side, cham_store_t storev, int M, int N, int K, int L, const CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, CHAMELEON_Complex64_t * W, int LDW) = NULL;
int (*libCORE_zparfb)(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zpemv)(cham_trans_t trans, cham_store_t storev, int M, int N, int L, CHAMELEON_Complex64_t ALPHA, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * X, int INCX, CHAMELEON_Complex64_t BETA, CHAMELEON_Complex64_t * Y, int INCY, CHAMELEON_Complex64_t * WORK) = NULL;
#if defined(PRECISION_z) || defined(PRECISION_c)
void (*libCORE_zplghe)(double bump, int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
#endif
void (*libCORE_zplgsy)(CHAMELEON_Complex64_t bump, int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
void (*libCORE_zplrnt)(int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) = NULL;
int (*libCORE_zplssq)(cham_store_t storev, int M, int N, double *sclssqin, double *sclssqout) = NULL;
int (*libCORE_zplssq2)(int N, double *sclssq) = NULL;
void (*libCORE_zpotrf)(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA, int *INFO) = NULL;
void (*libCORE_zshift)(int s, int m, int n, int L, CHAMELEON_Complex64_t * A) = NULL;
void (*libCORE_zshiftw)(int s, int cl, int m, int n, int L, CHAMELEON_Complex64_t * A, CHAMELEON_Complex64_t * W) = NULL;
int (*libCORE_zssssm)(int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * L1, int LDL1, const CHAMELEON_Complex64_t * L2, int LDL2, const int *IPIV) = NULL;
void (*libCORE_zsymm)(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
void (*libCORE_zsyrk)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
void (*libCORE_zsyr2k)(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) = NULL;
int (*libCORE_zsyssq)(cham_store_t storev, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) = NULL;
int (*libCORE_zsytf2_nopiv)(cham_uplo_t uplo, int n, CHAMELEON_Complex64_t * A, int lda) = NULL;
int (*libCORE_ztradd)(cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * B, int LDB) = NULL;
void (*libCORE_ztrasm)(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int lda, double *work) = NULL;
void (*libCORE_ztrmm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
void (*libCORE_ztrsm)(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) = NULL;
int (*libCORE_ztrssq)(cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *scale, double *sumsq) = NULL;
void (*libCORE_ztrtri)(cham_uplo_t uplo, cham_diag_t diag, int N, CHAMELEON_Complex64_t * A, int LDA, int *info) = NULL;
int (*libCORE_ztslqt)(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_ztsmlq)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_ztsmlq_hetra1)(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAMELEON_Complex64_t * A1, int lda1, CHAMELEON_Complex64_t * A2, int lda2, const CHAMELEON_Complex64_t * V, int ldv, const CHAMELEON_Complex64_t * T, int ldt, CHAMELEON_Complex64_t * WORK, int ldwork) = NULL;
int (*libCORE_ztsmqr)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_ztsmqr_hetra1)(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAMELEON_Complex64_t * A1, int lda1, CHAMELEON_Complex64_t * A2, int lda2, const CHAMELEON_Complex64_t * V, int ldv, const CHAMELEON_Complex64_t * T, int ldt, CHAMELEON_Complex64_t * WORK, int ldwork) = NULL;
int (*libCORE_ztsqrt)(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_ztstrf)(int M, int N, int IB, int NB, CHAMELEON_Complex64_t * U, int LDU, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * L, int LDL, int *IPIV, CHAMELEON_Complex64_t * WORK, int LDWORK, int *INFO) = NULL;
int (*libCORE_ztplqt)(int M, int N, int L, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_ztpmlqt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_ztpmqrt)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_ztpqrt)(int M, int N, int L, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_zttmqr)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zttqrt)(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_zttmlq)(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zttlqt)(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) = NULL;
int (*libCORE_zunmlq)(cham_side_t side, cham_trans_t trans, int M, int N, int IB, int K, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zunmqr)(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) = NULL;
int (*libCORE_zgesum)(cham_store_t storev, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * sum) = NULL;
int (*libCORE_zcesca)(int center, int scale, cham_store_t axis, int M, int N, int Mt, int Nt, const CHAMELEON_Complex64_t * Gi, int LDGI, const CHAMELEON_Complex64_t * Gj, int LDGJ, const CHAMELEON_Complex64_t * G, const double *Di, int LDDI, const double *Dj, int LDDJ, CHAMELEON_Complex64_t * A, int LDA) = NULL;
int (*libCORE_zgram)(cham_uplo_t uplo, int M, int N, int Mt, int Nt, const double *Di, int LDDI, const double *Dj, int LDDJ, const double *D, double *A, int LDA) = NULL;
void (*libCORE_zprint)(FILE * file, const char *header, cham_uplo_t uplo, int m, int n, const CHAMELEON_Complex64_t * A, int lda) = NULL;



/* Wrapper functions */


void CORE_dlag2z(cham_uplo_t uplo, int M, int N, const double *A, int lda, CHAMELEON_Complex64_t * B, int ldb) {
    FUNCTION_ENTRY;
    libCORE_dlag2z(uplo, M, N, A, lda, B, ldb);
    FUNCTION_EXIT;

}

void CORE_dzasum(cham_store_t storev, cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int lda, double *work) {
    FUNCTION_ENTRY;
    libCORE_dzasum(storev, uplo, M, N, A, lda, work);
    FUNCTION_EXIT;

}

int CORE_zaxpy(int M, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int incA, CHAMELEON_Complex64_t * B, int incB) {
    FUNCTION_ENTRY;
    int ret = libCORE_zaxpy(M, alpha, A, incA, B, incB);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgeadd(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgeadd(trans, M, N, alpha, A, LDA, beta, B, LDB);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlascal(cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t * A, int lda) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlascal(uplo, m, n, alpha, A, lda);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgelqt(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgelqt(M, N, IB, A, LDA, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zgemm(cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zgemm(transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    FUNCTION_EXIT;

}

void CORE_zgemv(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * x, int incX, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * y, int incY) {
    FUNCTION_ENTRY;
    libCORE_zgemv(trans, M, N, alpha, A, LDA, x, incX, beta, y, incY);
    FUNCTION_EXIT;

}

int CORE_zgeqrt(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgeqrt(M, N, IB, A, LDA, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgesplit(cham_side_t side, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgesplit(side, diag, M, N, A, LDA, B, LDB);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgessm(int M, int N, int K, int IB, const int *IPIV, const CHAMELEON_Complex64_t * L, int LDL, CHAMELEON_Complex64_t * A, int LDA) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgessm(M, N, K, IB, IPIV, L, LDL, A, LDA);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgessq(cham_store_t storev, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgessq(storev, M, N, A, LDA, sclssq);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetf2_nopiv(int M, int N, CHAMELEON_Complex64_t * A, int LDA) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetf2_nopiv(M, N, A, LDA);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf(int M, int N, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf(M, N, A, LDA, IPIV, INFO);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_incpiv(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_incpiv(M, N, IB, A, LDA, IPIV, INFO);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_nopiv(int M, int N, int IB, CHAMELEON_Complex64_t * A, int LDA, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_nopiv(M, N, IB, A, LDA, INFO);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_panel_diag(int m, int n, int h, int m0, int ib, CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * U, int ldu, int *IPIV, CHAM_pivot_t * nextpiv, const CHAM_pivot_t * prevpiv) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_panel_diag(m, n, h, m0, ib, A, lda, U, ldu, IPIV, nextpiv, prevpiv);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_panel_offdiag(int m, int n, int h, int m0, int ib, CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * U, int ldu, CHAM_pivot_t * nextpiv, const CHAM_pivot_t * prevpiv) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_panel_offdiag(m, n, h, m0, ib, A, lda, U, ldu, nextpiv, prevpiv);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_reclap(int M, int N, CHAMELEON_Complex64_t * A, int LDA, int *IPIV, int *info) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_reclap(M, N, A, LDA, IPIV, info);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgetrf_rectil(const CHAM_desc_t A, int *IPIV, int *info) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgetrf_rectil(A, IPIV, info);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zgetrip(int m, int n, CHAMELEON_Complex64_t * A, CHAMELEON_Complex64_t * work) {
    FUNCTION_ENTRY;
    libCORE_zgetrip(m, n, A, work);
    FUNCTION_EXIT;

}

void CORE_zhe2ge(cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    libCORE_zhe2ge(uplo, M, N, A, LDA, B, LDB);
    FUNCTION_EXIT;

}

void CORE_zhegst(int itype, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, int *INFO) {
    FUNCTION_ENTRY;
    libCORE_zhegst(itype, uplo, N, A, LDA, B, LDB, INFO);
    FUNCTION_EXIT;

}
#if defined(PRECISION_z) || defined(PRECISION_c)

void CORE_zhemm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zhemm(side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    FUNCTION_EXIT;

}

void CORE_zherk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, const CHAMELEON_Complex64_t * A, int LDA, double beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zherk(uplo, trans, N, K, alpha, A, LDA, beta, C, LDC);
    FUNCTION_EXIT;

}

void CORE_zher2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, double beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zher2k(uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    FUNCTION_EXIT;

}

int CORE_zhessq(cham_store_t storev, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) {
    FUNCTION_ENTRY;
    int ret = libCORE_zhessq(storev, uplo, N, A, LDA, sclssq);
    FUNCTION_EXIT;
    return ret;
}
#endif

int CORE_zherfb(cham_uplo_t uplo, int N, int K, int IB, int NB, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zherfb(uplo, N, K, IB, NB, A, LDA, T, LDT, C, LDC, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zlacpy(cham_uplo_t uplo, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    libCORE_zlacpy(uplo, M, N, A, LDA, B, LDB);
    FUNCTION_EXIT;

}

void CORE_zlange(cham_normtype_t norm, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) {
    FUNCTION_ENTRY;
    libCORE_zlange(norm, M, N, A, LDA, work, normA);
    FUNCTION_EXIT;

}
#if defined(PRECISION_z) || defined(PRECISION_c)

void CORE_zlanhe(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) {
    FUNCTION_ENTRY;
    libCORE_zlanhe(norm, uplo, N, A, LDA, work, normA);
    FUNCTION_EXIT;

}
#endif

void CORE_zlansy(cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) {
    FUNCTION_ENTRY;
    libCORE_zlansy(norm, uplo, N, A, LDA, work, normA);
    FUNCTION_EXIT;

}

void CORE_zlantr(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *work, double *normA) {
    FUNCTION_ENTRY;
    libCORE_zlantr(norm, uplo, diag, M, N, A, LDA, work, normA);
    FUNCTION_EXIT;

}

int CORE_zlarfb_gemm(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlarfb_gemm(side, trans, direct, storev, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlarfx2(cham_side_t side, int N, CHAMELEON_Complex64_t V, CHAMELEON_Complex64_t TAU, CHAMELEON_Complex64_t * C1, int LDC1, CHAMELEON_Complex64_t * C2, int LDC2) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlarfx2(side, N, V, TAU, C1, LDC1, C2, LDC2);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlarfx2c(cham_uplo_t uplo, CHAMELEON_Complex64_t V, CHAMELEON_Complex64_t TAU, CHAMELEON_Complex64_t * C1, CHAMELEON_Complex64_t * C2, CHAMELEON_Complex64_t * C3) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlarfx2c(uplo, V, TAU, C1, C2, C3);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlarfx2ce(cham_uplo_t uplo, CHAMELEON_Complex64_t * V, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * C1, CHAMELEON_Complex64_t * C2, CHAMELEON_Complex64_t * C3) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlarfx2ce(uplo, V, TAU, C1, C2, C3);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zlarfy(int N, CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * V, const CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    libCORE_zlarfy(N, A, LDA, V, TAU, WORK);
    FUNCTION_EXIT;

}

void CORE_zlaset(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * tileA, int ldtilea) {
    FUNCTION_ENTRY;
    libCORE_zlaset(uplo, n1, n2, alpha, beta, tileA, ldtilea);
    FUNCTION_EXIT;

}

void CORE_zlaset2(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t * tileA, int ldtilea) {
    FUNCTION_ENTRY;
    libCORE_zlaset2(uplo, n1, n2, alpha, tileA, ldtilea);
    FUNCTION_EXIT;

}

void CORE_zlaswp(int N, CHAMELEON_Complex64_t * A, int LDA, int I1, int I2, const int *IPIV, int INC) {
    FUNCTION_ENTRY;
    libCORE_zlaswp(N, A, LDA, I1, I2, IPIV, INC);
    FUNCTION_EXIT;

}

int CORE_zlaswp_get(cham_side_t side, int m0, int m, int n, int k, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb, const int *perm) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlaswp_get(side, m0, m, n, k, A, lda, B, ldb, perm);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlaswp_set(cham_side_t side, int m0, int m, int n, int k, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb, const int *invp) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlaswp_set(side, m0, m, n, k, A, lda, B, ldb, invp);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlaswp_ontile(CHAM_desc_t descA, int i1, int i2, const int *ipiv, int inc) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlatro(cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlatro(uplo, trans, M, N, A, LDA, B, LDB);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zlatm1(int MODE, double COND, int IRSIGN, cham_dist_t DIST, unsigned long long int seed, CHAMELEON_Complex64_t * D, int N) {
    FUNCTION_ENTRY;
    int ret = libCORE_zlatm1(MODE, COND, IRSIGN, DIST, seed, D, N);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zlauum(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA) {
    FUNCTION_ENTRY;
    libCORE_zlauum(uplo, N, A, LDA);
    FUNCTION_EXIT;

}

int CORE_zpamm(int op, cham_side_t side, cham_store_t storev, int M, int N, int K, int L, const CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, CHAMELEON_Complex64_t * W, int LDW) {
    FUNCTION_ENTRY;
    int ret = libCORE_zpamm(op, side, storev, M, N, K, L, A1, LDA1, A2, LDA2, V, LDV, W, LDW);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zparfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zparfb(side, trans, direct, storev, M1, N1, M2, N2, K, L, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zpemv(cham_trans_t trans, cham_store_t storev, int M, int N, int L, CHAMELEON_Complex64_t ALPHA, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * X, int INCX, CHAMELEON_Complex64_t BETA, CHAMELEON_Complex64_t * Y, int INCY, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zpemv(trans, storev, M, N, L, ALPHA, A, LDA, X, INCX, BETA, Y, INCY, WORK);
    FUNCTION_EXIT;
    return ret;
}
#if defined(PRECISION_z) || defined(PRECISION_c)

void CORE_zplghe(double bump, int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libCORE_zplghe(bump, m, n, A, lda, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}
#endif

void CORE_zplgsy(CHAMELEON_Complex64_t bump, int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libCORE_zplgsy(bump, m, n, A, lda, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}

void CORE_zplrnt(int m, int n, CHAMELEON_Complex64_t * A, int lda, int bigM, int m0, int n0, unsigned long long int seed) {
    FUNCTION_ENTRY;
    libCORE_zplrnt(m, n, A, lda, bigM, m0, n0, seed);
    FUNCTION_EXIT;

}

int CORE_zplssq(cham_store_t storev, int M, int N, double *sclssqin, double *sclssqout) {
    FUNCTION_ENTRY;
    int ret = libCORE_zplssq(storev, M, N, sclssqin, sclssqout);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zplssq2(int N, double *sclssq) {
    FUNCTION_ENTRY;
    int ret = libCORE_zplssq2(N, sclssq);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zpotrf(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t * A, int LDA, int *INFO) {
    FUNCTION_ENTRY;
    libCORE_zpotrf(uplo, N, A, LDA, INFO);
    FUNCTION_EXIT;

}

void CORE_zshift(int s, int m, int n, int L, CHAMELEON_Complex64_t * A) {
    FUNCTION_ENTRY;
    libCORE_zshift(s, m, n, L, A);
    FUNCTION_EXIT;

}

void CORE_zshiftw(int s, int cl, int m, int n, int L, CHAMELEON_Complex64_t * A, CHAMELEON_Complex64_t * W) {
    FUNCTION_ENTRY;
    libCORE_zshiftw(s, cl, m, n, L, A, W);
    FUNCTION_EXIT;

}

int CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * L1, int LDL1, const CHAMELEON_Complex64_t * L2, int LDL2, const int *IPIV) {
    FUNCTION_ENTRY;
    int ret = libCORE_zssssm(M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, L1, LDL1, L2, LDL2, IPIV);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zsymm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zsymm(side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    FUNCTION_EXIT;

}

void CORE_zsyrk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zsyrk(uplo, trans, N, K, alpha, A, LDA, beta, C, LDC);
    FUNCTION_EXIT;

}

void CORE_zsyr2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, const CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * C, int LDC) {
    FUNCTION_ENTRY;
    libCORE_zsyr2k(uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    FUNCTION_EXIT;

}

int CORE_zsyssq(cham_store_t storev, cham_uplo_t uplo, int N, const CHAMELEON_Complex64_t * A, int LDA, double *sclssq) {
    FUNCTION_ENTRY;
    int ret = libCORE_zsyssq(storev, uplo, N, A, LDA, sclssq);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zsytf2_nopiv(cham_uplo_t uplo, int n, CHAMELEON_Complex64_t * A, int lda) {
    FUNCTION_ENTRY;
    int ret = libCORE_zsytf2_nopiv(uplo, n, A, lda);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztradd(cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    FUNCTION_EXIT;
    return ret;
}

void CORE_ztrasm(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int lda, double *work) {
    FUNCTION_ENTRY;
    libCORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
    FUNCTION_EXIT;

}

void CORE_ztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    libCORE_ztrmm(side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    FUNCTION_EXIT;

}

void CORE_ztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB) {
    FUNCTION_ENTRY;
    libCORE_ztrsm(side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    FUNCTION_EXIT;

}

int CORE_ztrssq(cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, double *scale, double *sumsq) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztrssq(uplo, diag, M, N, A, LDA, scale, sumsq);
    FUNCTION_EXIT;
    return ret;
}

void CORE_ztrtri(cham_uplo_t uplo, cham_diag_t diag, int N, CHAMELEON_Complex64_t * A, int LDA, int *info) {
    FUNCTION_ENTRY;
    libCORE_ztrtri(uplo, diag, N, A, LDA, info);
    FUNCTION_EXIT;

}

int CORE_ztslqt(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztslqt(M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztsmlq(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztsmlq(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztsmlq_hetra1(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAMELEON_Complex64_t * A1, int lda1, CHAMELEON_Complex64_t * A2, int lda2, const CHAMELEON_Complex64_t * V, int ldv, const CHAMELEON_Complex64_t * T, int ldt, CHAMELEON_Complex64_t * WORK, int ldwork) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztsmlq_hetra1(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztsmqr(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztsmqr(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztsmqr_hetra1(cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAMELEON_Complex64_t * A1, int lda1, CHAMELEON_Complex64_t * A2, int lda2, const CHAMELEON_Complex64_t * V, int ldv, const CHAMELEON_Complex64_t * T, int ldt, CHAMELEON_Complex64_t * WORK, int ldwork) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztsmqr_hetra1(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztsqrt(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztsqrt(M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztstrf(int M, int N, int IB, int NB, CHAMELEON_Complex64_t * U, int LDU, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * L, int LDL, int *IPIV, CHAMELEON_Complex64_t * WORK, int LDWORK, int *INFO) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztstrf(M, N, IB, NB, U, LDU, A, LDA, L, LDL, IPIV, WORK, LDWORK, INFO);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztplqt(int M, int N, int L, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztplqt(M, N, L, IB, A, LDA, B, LDB, T, LDT, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztpmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztpmlqt(side, trans, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztpmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztpmqrt(side, trans, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_ztpqrt(int M, int N, int L, int IB, CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * B, int LDB, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_ztpqrt(M, N, L, IB, A, LDA, B, LDB, T, LDT, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zttmqr(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zttmqr(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zttqrt(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zttqrt(M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zttmlq(cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zttmlq(side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zttlqt(int M, int N, int IB, CHAMELEON_Complex64_t * A1, int LDA1, CHAMELEON_Complex64_t * A2, int LDA2, CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * TAU, CHAMELEON_Complex64_t * WORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zttlqt(M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zunmlq(cham_side_t side, cham_trans_t trans, int M, int N, int IB, int K, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zunmlq(side, trans, M, N, IB, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zunmqr(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAMELEON_Complex64_t * V, int LDV, const CHAMELEON_Complex64_t * T, int LDT, CHAMELEON_Complex64_t * C, int LDC, CHAMELEON_Complex64_t * WORK, int LDWORK) {
    FUNCTION_ENTRY;
    int ret = libCORE_zunmqr(side, trans, M, N, K, IB, V, LDV, T, LDT, C, LDC, WORK, LDWORK);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgesum(cham_store_t storev, int M, int N, const CHAMELEON_Complex64_t * A, int LDA, CHAMELEON_Complex64_t * sum) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgesum(storev, M, N, A, LDA, sum);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zcesca(int center, int scale, cham_store_t axis, int M, int N, int Mt, int Nt, const CHAMELEON_Complex64_t * Gi, int LDGI, const CHAMELEON_Complex64_t * Gj, int LDGJ, const CHAMELEON_Complex64_t * G, const double *Di, int LDDI, const double *Dj, int LDDJ, CHAMELEON_Complex64_t * A, int LDA) {
    FUNCTION_ENTRY;
    int ret = libCORE_zcesca(center, scale, axis, M, N, Mt, Nt, Gi, LDGI, Gj, LDGJ, G, Di, LDDI, Dj, LDDJ, A, LDA);
    FUNCTION_EXIT;
    return ret;
}

int CORE_zgram(cham_uplo_t uplo, int M, int N, int Mt, int Nt, const double *Di, int LDDI, const double *Dj, int LDDJ, const double *D, double *A, int LDA) {
    FUNCTION_ENTRY;
    int ret = libCORE_zgram(uplo, M, N, Mt, Nt, Di, LDDI, Dj, LDDJ, D, A, LDA);
    FUNCTION_EXIT;
    return ret;
}

void CORE_zprint(FILE * file, const char *header, cham_uplo_t uplo, int m, int n, const CHAMELEON_Complex64_t * A, int lda) {
    FUNCTION_ENTRY;
    libCORE_zprint(file, header, uplo, m, n, A, lda);
    FUNCTION_EXIT;

}


PPTRACE_START_INTERCEPT_FUNCTIONS(coreblas_core_z)
INTERCEPT3("CORE_dlag2z", libCORE_dlag2z)
INTERCEPT3("CORE_dzasum", libCORE_dzasum)
INTERCEPT3("CORE_zaxpy", libCORE_zaxpy)
INTERCEPT3("CORE_zgeadd", libCORE_zgeadd)
INTERCEPT3("CORE_zlascal", libCORE_zlascal)
INTERCEPT3("CORE_zgelqt", libCORE_zgelqt)
INTERCEPT3("CORE_zgemm", libCORE_zgemm)
INTERCEPT3("CORE_zgemv", libCORE_zgemv)
INTERCEPT3("CORE_zgeqrt", libCORE_zgeqrt)
INTERCEPT3("CORE_zgesplit", libCORE_zgesplit)
INTERCEPT3("CORE_zgessm", libCORE_zgessm)
INTERCEPT3("CORE_zgessq", libCORE_zgessq)
INTERCEPT3("CORE_zgetf2_nopiv", libCORE_zgetf2_nopiv)
INTERCEPT3("CORE_zgetrf", libCORE_zgetrf)
INTERCEPT3("CORE_zgetrf_incpiv", libCORE_zgetrf_incpiv)
INTERCEPT3("CORE_zgetrf_nopiv", libCORE_zgetrf_nopiv)
INTERCEPT3("CORE_zgetrf_panel_diag", libCORE_zgetrf_panel_diag)
INTERCEPT3("CORE_zgetrf_panel_offdiag", libCORE_zgetrf_panel_offdiag)
INTERCEPT3("CORE_zgetrf_reclap", libCORE_zgetrf_reclap)
INTERCEPT3("CORE_zgetrf_rectil", libCORE_zgetrf_rectil)
INTERCEPT3("CORE_zgetrip", libCORE_zgetrip)
INTERCEPT3("CORE_zhe2ge", libCORE_zhe2ge)
INTERCEPT3("CORE_zhegst", libCORE_zhegst)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("CORE_zhemm", libCORE_zhemm)
INTERCEPT3("CORE_zherk", libCORE_zherk)
INTERCEPT3("CORE_zher2k", libCORE_zher2k)
INTERCEPT3("CORE_zhessq", libCORE_zhessq)
#endif
INTERCEPT3("CORE_zherfb", libCORE_zherfb)
INTERCEPT3("CORE_zlacpy", libCORE_zlacpy)
INTERCEPT3("CORE_zlange", libCORE_zlange)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("CORE_zlanhe", libCORE_zlanhe)
#endif
INTERCEPT3("CORE_zlansy", libCORE_zlansy)
INTERCEPT3("CORE_zlantr", libCORE_zlantr)
INTERCEPT3("CORE_zlarfb_gemm", libCORE_zlarfb_gemm)
INTERCEPT3("CORE_zlarfx2", libCORE_zlarfx2)
INTERCEPT3("CORE_zlarfx2c", libCORE_zlarfx2c)
INTERCEPT3("CORE_zlarfx2ce", libCORE_zlarfx2ce)
INTERCEPT3("CORE_zlarfy", libCORE_zlarfy)
INTERCEPT3("CORE_zlaset", libCORE_zlaset)
INTERCEPT3("CORE_zlaset2", libCORE_zlaset2)
INTERCEPT3("CORE_zlaswp", libCORE_zlaswp)
INTERCEPT3("CORE_zlaswp_get", libCORE_zlaswp_get)
INTERCEPT3("CORE_zlaswp_set", libCORE_zlaswp_set)
INTERCEPT3("CORE_zlaswp_ontile", libCORE_zlaswp_ontile)
INTERCEPT3("CORE_zlatro", libCORE_zlatro)
INTERCEPT3("CORE_zlatm1", libCORE_zlatm1)
INTERCEPT3("CORE_zlauum", libCORE_zlauum)
INTERCEPT3("CORE_zpamm", libCORE_zpamm)
INTERCEPT3("CORE_zparfb", libCORE_zparfb)
INTERCEPT3("CORE_zpemv", libCORE_zpemv)
#if defined(PRECISION_z) || defined(PRECISION_c)
INTERCEPT3("CORE_zplghe", libCORE_zplghe)
#endif
INTERCEPT3("CORE_zplgsy", libCORE_zplgsy)
INTERCEPT3("CORE_zplrnt", libCORE_zplrnt)
INTERCEPT3("CORE_zplssq", libCORE_zplssq)
INTERCEPT3("CORE_zplssq2", libCORE_zplssq2)
INTERCEPT3("CORE_zpotrf", libCORE_zpotrf)
INTERCEPT3("CORE_zshift", libCORE_zshift)
INTERCEPT3("CORE_zshiftw", libCORE_zshiftw)
INTERCEPT3("CORE_zssssm", libCORE_zssssm)
INTERCEPT3("CORE_zsymm", libCORE_zsymm)
INTERCEPT3("CORE_zsyrk", libCORE_zsyrk)
INTERCEPT3("CORE_zsyr2k", libCORE_zsyr2k)
INTERCEPT3("CORE_zsyssq", libCORE_zsyssq)
INTERCEPT3("CORE_zsytf2_nopiv", libCORE_zsytf2_nopiv)
INTERCEPT3("CORE_ztradd", libCORE_ztradd)
INTERCEPT3("CORE_ztrasm", libCORE_ztrasm)
INTERCEPT3("CORE_ztrmm", libCORE_ztrmm)
INTERCEPT3("CORE_ztrsm", libCORE_ztrsm)
INTERCEPT3("CORE_ztrssq", libCORE_ztrssq)
INTERCEPT3("CORE_ztrtri", libCORE_ztrtri)
INTERCEPT3("CORE_ztslqt", libCORE_ztslqt)
INTERCEPT3("CORE_ztsmlq", libCORE_ztsmlq)
INTERCEPT3("CORE_ztsmlq_hetra1", libCORE_ztsmlq_hetra1)
INTERCEPT3("CORE_ztsmqr", libCORE_ztsmqr)
INTERCEPT3("CORE_ztsmqr_hetra1", libCORE_ztsmqr_hetra1)
INTERCEPT3("CORE_ztsqrt", libCORE_ztsqrt)
INTERCEPT3("CORE_ztstrf", libCORE_ztstrf)
INTERCEPT3("CORE_ztplqt", libCORE_ztplqt)
INTERCEPT3("CORE_ztpmlqt", libCORE_ztpmlqt)
INTERCEPT3("CORE_ztpmqrt", libCORE_ztpmqrt)
INTERCEPT3("CORE_ztpqrt", libCORE_ztpqrt)
INTERCEPT3("CORE_zttmqr", libCORE_zttmqr)
INTERCEPT3("CORE_zttqrt", libCORE_zttqrt)
INTERCEPT3("CORE_zttmlq", libCORE_zttmlq)
INTERCEPT3("CORE_zttlqt", libCORE_zttlqt)
INTERCEPT3("CORE_zunmlq", libCORE_zunmlq)
INTERCEPT3("CORE_zunmqr", libCORE_zunmqr)
INTERCEPT3("CORE_zgesum", libCORE_zgesum)
INTERCEPT3("CORE_zcesca", libCORE_zcesca)
INTERCEPT3("CORE_zgram", libCORE_zgram)
INTERCEPT3("CORE_zprint", libCORE_zprint)

PPTRACE_END_INTERCEPT_FUNCTIONS(coreblas_core_z)
