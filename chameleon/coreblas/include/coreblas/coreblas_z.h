/**
 *
 * @file coreblas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @author Brieuc Nicolas
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _coreblas_z_h_
#define _coreblas_z_h_

/**
 *  Declarations of serial kernels - alphabetical order
 */
void CORE_dlag2z( cham_uplo_t uplo, int M, int N,
                  const double          *A, int lda,
                  CHAMELEON_Complex64_t *B, int ldb );
void CORE_dzasum(cham_store_t storev, cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int lda, double *work);
int CORE_zaxpy(int M, CHAMELEON_Complex64_t alpha,
               const CHAMELEON_Complex64_t *A, int incA,
                     CHAMELEON_Complex64_t *B, int incB);
int CORE_zgeadd(cham_trans_t trans, int M, int N,
                      CHAMELEON_Complex64_t alpha,
                const CHAMELEON_Complex64_t *A, int LDA,
                      CHAMELEON_Complex64_t beta,
                      CHAMELEON_Complex64_t *B, int LDB);
int CORE_zlascal( cham_uplo_t uplo, int m, int n,
              CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int lda );
int  CORE_zgelqt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU,
                 CHAMELEON_Complex64_t *WORK);
void CORE_zgemm(cham_trans_t transA, cham_trans_t transB,
                int M, int N, int K,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                          const CHAMELEON_Complex64_t *B, int LDB,
                CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *C, int LDC);
void CORE_zgemv(cham_trans_t trans, int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                             const CHAMELEON_Complex64_t *x, int incX,
                CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *y, int incY);
int  CORE_zgeqrt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK);
int CORE_zgesplit(cham_side_t side, cham_diag_t diag,
                  int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t *B, int LDB);
int  CORE_zgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const CHAMELEON_Complex64_t *L, int LDL,
                 CHAMELEON_Complex64_t *A, int LDA);
int  CORE_zgessq(cham_store_t storev, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq);
int CORE_zgetf2_nopiv(int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA);
int  CORE_zgetrf(int M, int N,
                 CHAMELEON_Complex64_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_zgetrf_incpiv(int M, int N, int IB,
                        CHAMELEON_Complex64_t *A, int LDA,
                        int *IPIV, int *INFO);
int CORE_zgetrf_nopiv(int M, int N, int IB,
                      CHAMELEON_Complex64_t *A, int LDA,
                      int *INFO);
int CORE_zgetrf_panel_diag( int m, int n, int h, int m0, int ib,
                            CHAMELEON_Complex64_t *A, int lda,
                            CHAMELEON_Complex64_t *U, int ldu,
                            int *IPIV,
                            CHAM_pivot_t *nextpiv, const CHAM_pivot_t *prevpiv);
int CORE_zgetrf_panel_offdiag( int m, int n, int h, int m0, int ib,
                               CHAMELEON_Complex64_t *A, int lda,
                               CHAMELEON_Complex64_t *U, int ldu,
                               CHAM_pivot_t *nextpiv, const CHAM_pivot_t *prevpiv );
int  CORE_zgetrf_reclap(int M, int N,
                        CHAMELEON_Complex64_t *A, int LDA,
                        int *IPIV, int *info);
int  CORE_zgetrf_rectil(const CHAM_desc_t A, int *IPIV, int *info);
void CORE_zgetrip(int m, int n, CHAMELEON_Complex64_t *A,
                  CHAMELEON_Complex64_t *work);
void CORE_zhe2ge(cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB);
void CORE_zhegst(int itype, cham_uplo_t uplo, int N,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB, int *INFO);
#if defined(PRECISION_z) || defined(PRECISION_c)
void CORE_zhemm(cham_side_t side, cham_uplo_t uplo,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                          const CHAMELEON_Complex64_t *B, int LDB,
                CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *C, int LDC);
void CORE_zherk(cham_uplo_t uplo, cham_trans_t trans,
                int N, int K,
                double alpha, const CHAMELEON_Complex64_t *A, int LDA,
                double beta,        CHAMELEON_Complex64_t *C, int LDC);
void CORE_zher2k(cham_uplo_t uplo, cham_trans_t trans,
                 int N, int K,
                 CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                           const CHAMELEON_Complex64_t *B, int LDB,
                 double beta,                    CHAMELEON_Complex64_t *C, int LDC);
int  CORE_zhessq(cham_store_t storev, cham_uplo_t uplo, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq);
#endif
int  CORE_zherfb(cham_uplo_t uplo, int N, int K, int IB, int NB,
                 const CHAMELEON_Complex64_t *A,    int LDA,
                 const CHAMELEON_Complex64_t *T,    int LDT,
                       CHAMELEON_Complex64_t *C,    int LDC,
                       CHAMELEON_Complex64_t *WORK, int LDWORK);
void CORE_zlacpy(cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                       CHAMELEON_Complex64_t *B, int LDB);
void CORE_zlange(cham_normtype_t norm, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *work, double *normA);
#if defined(PRECISION_z) || defined(PRECISION_c)
void CORE_zlanhe(cham_normtype_t norm, cham_uplo_t uplo, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *work, double *normA);
#endif
void CORE_zlansy(cham_normtype_t norm, cham_uplo_t uplo, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlantr(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *work, double *normA);
int CORE_zlarfb_gemm(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev,
                     int M, int N, int K,
                     const CHAMELEON_Complex64_t *V, int LDV,
                     const CHAMELEON_Complex64_t *T, int LDT,
                           CHAMELEON_Complex64_t *C, int LDC,
                           CHAMELEON_Complex64_t *WORK, int LDWORK);
int CORE_zlarfx2(cham_side_t side, int N,
                 CHAMELEON_Complex64_t V,
                 CHAMELEON_Complex64_t TAU,
                 CHAMELEON_Complex64_t *C1, int LDC1,
                 CHAMELEON_Complex64_t *C2, int LDC2);
int CORE_zlarfx2c(cham_uplo_t uplo,
                  CHAMELEON_Complex64_t V,
                  CHAMELEON_Complex64_t TAU,
                  CHAMELEON_Complex64_t *C1,
                  CHAMELEON_Complex64_t *C2,
                  CHAMELEON_Complex64_t *C3);
int CORE_zlarfx2ce(cham_uplo_t uplo,
                   CHAMELEON_Complex64_t *V,
                   CHAMELEON_Complex64_t *TAU,
                   CHAMELEON_Complex64_t *C1,
                   CHAMELEON_Complex64_t *C2,
                   CHAMELEON_Complex64_t *C3);
void CORE_zlarfy(int N,
                 CHAMELEON_Complex64_t *A, int LDA,
                 const CHAMELEON_Complex64_t *V,
                 const CHAMELEON_Complex64_t *TAU,
                 CHAMELEON_Complex64_t *WORK);
void CORE_zlaset(cham_uplo_t uplo, int n1, int n2,
                 CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                 CHAMELEON_Complex64_t *tileA, int ldtilea);
void CORE_zlaset2(cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                  CHAMELEON_Complex64_t *tileA, int ldtilea);
void CORE_zlaswp(int N, CHAMELEON_Complex64_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int CORE_zlaswp_get( cham_side_t side, int m0, int m, int n, int k,
                     const CHAMELEON_Complex64_t *A, int lda,
                     CHAMELEON_Complex64_t *B, int ldb,
                     const int *perm );
int CORE_zlaswp_get_idx( cham_side_t side, int m0, int m, int n, int k,
                         const CHAMELEON_Complex64_t *A, int lda,
                         CHAMELEON_Complex64_t *B, int ldb,
                         int perm_m, int perm_mt, const int *perm_idx );
int CORE_zlaswp_set( cham_side_t side, int m0, int m, int n, int k,
                     const CHAMELEON_Complex64_t *A, int lda,
                     CHAMELEON_Complex64_t *B, int ldb,
                     const int *invp );
int CORE_zlaswp_set_idx( cham_side_t side, int m0, int m, int n, int k,
                         const CHAMELEON_Complex64_t *A, int lda,
                         CHAMELEON_Complex64_t *B, int ldb,
                         int invp_m, int invp_mt, const int *invp_idx );
int  CORE_zlaswp_ontile( CHAM_desc_t descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlatro(cham_uplo_t uplo, cham_trans_t trans,
                 int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                       CHAMELEON_Complex64_t *B, int LDB);
int CORE_zlatm1( int MODE, double COND, int IRSIGN, cham_dist_t DIST,
                 unsigned long long int seed,
                 CHAMELEON_Complex64_t *D, int N );
void CORE_zlauum(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CORE_zpamm(int op, cham_side_t side, cham_store_t storev,
               int M, int N, int K, int L,
               const CHAMELEON_Complex64_t *A1, int LDA1,
                     CHAMELEON_Complex64_t *A2, int LDA2,
               const CHAMELEON_Complex64_t *V, int LDV,
                     CHAMELEON_Complex64_t *W, int LDW);
int  CORE_zparfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       CHAMELEON_Complex64_t *A1, int LDA1,
                       CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                       CHAMELEON_Complex64_t *WORK, int LDWORK);
int CORE_zpemv(cham_trans_t trans, cham_store_t storev,
               int M, int N, int L,
               CHAMELEON_Complex64_t ALPHA,
               const CHAMELEON_Complex64_t *A, int LDA,
               const CHAMELEON_Complex64_t *X, int INCX,
               CHAMELEON_Complex64_t BETA,
               CHAMELEON_Complex64_t *Y, int INCY,
               CHAMELEON_Complex64_t *WORK);
#if defined(PRECISION_z) || defined(PRECISION_c)
void CORE_zplghe(double bump, int m, int n, CHAMELEON_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
#endif
void CORE_zplgsy(CHAMELEON_Complex64_t bump, int m, int n, CHAMELEON_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplrnt(int m, int n, CHAMELEON_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int CORE_zplssq( cham_store_t storev, int M, int N,
                 double *sclssqin, double *sclssqout );
int CORE_zplssq2( int N, double *sclssq );
void CORE_zpotrf(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, int *INFO);
void CORE_zshift(int s, int m, int n, int L,
                 CHAMELEON_Complex64_t *A);
void CORE_zshiftw(int s, int cl, int m, int n, int L,
                  CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *W);
int  CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       CHAMELEON_Complex64_t *A1, int LDA1,
                       CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *L1, int LDL1,
                 const CHAMELEON_Complex64_t *L2, int LDL2,
                 const int *IPIV);
void CORE_zsymm(cham_side_t side, cham_uplo_t uplo,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                          const CHAMELEON_Complex64_t *B, int LDB,
                CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *C, int LDC);
void CORE_zsyrk(cham_uplo_t uplo, cham_trans_t trans,
                int N, int K,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *C, int LDC);
void CORE_zsyr2k(cham_uplo_t uplo, cham_trans_t trans,
                 int N, int K,
                 CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                           const CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *C, int LDC);
int  CORE_zsyssq(cham_store_t storev, cham_uplo_t uplo, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq);
int CORE_zsytf2_nopiv(cham_uplo_t uplo, int n, CHAMELEON_Complex64_t *A, int lda);
int CORE_ztradd(cham_uplo_t uplo, cham_trans_t trans, int M, int N,
                      CHAMELEON_Complex64_t alpha,
                const CHAMELEON_Complex64_t *A, int LDA,
                      CHAMELEON_Complex64_t beta,
                      CHAMELEON_Complex64_t *B, int LDB);
void CORE_ztrasm(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag,
                 int M, int N,
                 const CHAMELEON_Complex64_t *A, int lda, double *work);
void CORE_ztrmm(cham_side_t side, cham_uplo_t uplo,
                cham_trans_t transA, cham_diag_t diag,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                                CHAMELEON_Complex64_t *B, int LDB);
void CORE_ztrsm(cham_side_t side, cham_uplo_t uplo,
                cham_trans_t transA, cham_diag_t diag,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                                CHAMELEON_Complex64_t *B, int LDB);
int CORE_ztrssq(cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                const CHAMELEON_Complex64_t *A, int LDA,
                double *scale, double *sumsq);
void CORE_ztrtri(cham_uplo_t uplo, cham_diag_t diag, int N,
                 CHAMELEON_Complex64_t *A, int LDA, int *info);
int  CORE_ztslqt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK);
int  CORE_ztsmlq(cham_side_t side, cham_trans_t trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int CORE_ztsmlq_hetra1( cham_side_t side, cham_trans_t trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        CHAMELEON_Complex64_t *A1, int lda1,
                        CHAMELEON_Complex64_t *A2, int lda2,
                        const CHAMELEON_Complex64_t *V, int ldv,
                        const CHAMELEON_Complex64_t *T, int ldt,
                        CHAMELEON_Complex64_t *WORK, int ldwork);
int  CORE_ztsmqr(cham_side_t side, cham_trans_t trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int CORE_ztsmqr_hetra1( cham_side_t side, cham_trans_t trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        CHAMELEON_Complex64_t *A1, int lda1,
                        CHAMELEON_Complex64_t *A2, int lda2,
                        const CHAMELEON_Complex64_t *V, int ldv,
                        const CHAMELEON_Complex64_t *T, int ldt,
                        CHAMELEON_Complex64_t *WORK, int ldwork);
int  CORE_ztsqrt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK);
int  CORE_ztstrf(int M, int N, int IB, int NB,
                 CHAMELEON_Complex64_t *U, int LDU,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *L, int LDL,
                 int *IPIV, CHAMELEON_Complex64_t *WORK,
                 int LDWORK, int *INFO);
int CORE_ztplqt( int M, int N, int L, int IB,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK );
int CORE_ztpmlqt( cham_side_t side, cham_trans_t trans,
                  int M, int N, int K, int L, int IB,
                  const CHAMELEON_Complex64_t *V, int LDV,
                  const CHAMELEON_Complex64_t *T, int LDT,
                  CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t *B, int LDB,
                  CHAMELEON_Complex64_t *WORK );
int CORE_ztpmqrt( cham_side_t side, cham_trans_t trans,
                  int M, int N, int K, int L, int IB,
                  const CHAMELEON_Complex64_t *V, int LDV,
                  const CHAMELEON_Complex64_t *T, int LDT,
                  CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t *B, int LDB,
                  CHAMELEON_Complex64_t *WORK );
int CORE_ztpqrt( int M, int N, int L, int IB,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK );
int  CORE_zttmqr(cham_side_t side, cham_trans_t trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int  CORE_zttqrt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU,
                 CHAMELEON_Complex64_t *WORK);
int  CORE_zttmlq(cham_side_t side, cham_trans_t trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int  CORE_zttlqt(int M, int N, int IB,
                 CHAMELEON_Complex64_t *A1, int LDA1,
                 CHAMELEON_Complex64_t *A2, int LDA2,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *TAU,
                 CHAMELEON_Complex64_t *WORK);
int  CORE_zunmlq(cham_side_t side, cham_trans_t trans,
                 int M, int N, int IB, int K,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *C, int LDC,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(cham_side_t side, cham_trans_t trans,
                 int M, int N, int K, int IB,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *C, int LDC,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);

/**
 * Centered-Scale and Gram prototypes
 */
int CORE_zgesum( cham_store_t storev, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *sum );
int CORE_zcesca( int center, int scale,
                 cham_store_t axis,
                 int M, int N,
                 int Mt, int Nt,
                 const CHAMELEON_Complex64_t *Gi, int LDGI,
                 const CHAMELEON_Complex64_t *Gj, int LDGJ,
                 const CHAMELEON_Complex64_t *G,
                 const double *Di, int LDDI,
                 const double *Dj, int LDDJ,
                 CHAMELEON_Complex64_t *A, int LDA );
int CORE_zgram( cham_uplo_t uplo, int M, int N, int Mt, int Nt,
                const double *Di, int LDDI,
                const double *Dj, int LDDJ,
                const double *D,
                double *A, int LDA );

void CORE_zprint( FILE *file, const char *header,
                  cham_uplo_t uplo, int m, int n,
                  const CHAMELEON_Complex64_t *A, int lda );

#endif /* _coreblas_z_h_ */
