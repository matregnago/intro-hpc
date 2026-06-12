/**
 *
 * @file chameleon_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_complex64_t wrappers
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
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Ana Hourcau
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_z_h_
#define _chameleon_z_h_

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
int CHAMELEON_zgeadd(cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *B, int LDB);
//int CHAMELEON_zgebrd(int M, int N, CHAMELEON_Complex64_t *A, int LDA, double *D, double *E, CHAM_desc_t *descT);
//int CHAMELEON_zgecon(cham_normtype_t norm, int N, CHAMELEON_Complex64_t *A, int LDA, double anorm, double *rcond);
//int CHAMELEON_zpocon(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, double anorm, double *rcond);
int CHAMELEON_zgelqf(int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT);
int CHAMELEON_zgelqs(int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgels(cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
double CHAMELEON_zgenm2( double tol, int M, int N, CHAMELEON_Complex64_t *A, int LDA );
int CHAMELEON_zgemm(cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zgepdf_qdwh( int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *H, int LDH, gepdf_info_t *info );
int CHAMELEON_zgeqrf(int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT);
int CHAMELEON_zgeqrs(int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgesv(int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, int *IPIV, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgesv_incpiv(int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descL, int *IPIV, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgesv_nopiv(int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgesvd(cham_job_t jobu, cham_job_t jobvt, int M, int N, CHAMELEON_Complex64_t *A, int LDA, double *S, CHAM_desc_t *descT, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT);
//int CHAMELEON_zgetrf(int M, int N, CHAMELEON_Complex64_t *A, int LDA, int *IPIV);
int CHAMELEON_zgetrf_incpiv(int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descL, int *IPIV);
int CHAMELEON_zgetrf_nopiv(int M, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zgetrf( int M, int N, CHAMELEON_Complex64_t *A, int LDA, int *IPIV );
//int CHAMELEON_zgetri(int N, CHAMELEON_Complex64_t *A, int LDA, int *IPIV);
int CHAMELEON_zgetrs(cham_trans_t trans, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, int *IPIV, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgetrs_incpiv(cham_trans_t trans, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descL, int *IPIV, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgetrs_nopiv(cham_trans_t trans, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zhemm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zherk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, CHAMELEON_Complex64_t *A, int LDA, double beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zher2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, double beta, CHAMELEON_Complex64_t *C, int LDC);
//int CHAMELEON_zheev(cham_job_t jobz, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, double *W, CHAM_desc_t *descT, CHAMELEON_Complex64_t *Q, int LDQ);
int CHAMELEON_zheevd(cham_job_t jobz, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, double *W, CHAM_desc_t *descT);
//int CHAMELEON_zhegv(CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, double *W, CHAM_desc_t *descT, CHAMELEON_Complex64_t *Q, int LDQ);
//int CHAMELEON_zhegvd(CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, double *W, CHAM_desc_t *descT, CHAMELEON_Complex64_t *Q, int LDQ);
//int CHAMELEON_zhegst(CHAMELEON_enum itype, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zhetrd(cham_job_t jobz, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, double *D, double *E, CHAM_desc_t *descT, CHAMELEON_Complex64_t *Q, int LDQ);
int CHAMELEON_zlacpy(cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
double CHAMELEON_zlange(cham_normtype_t norm, int M, int N, CHAMELEON_Complex64_t *A, int LDA);
double CHAMELEON_zlanhe(cham_normtype_t norm, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
double CHAMELEON_zlansy(cham_normtype_t norm, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
double CHAMELEON_zlantr(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zlascal(cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zlaset(cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zlaswp( cham_side_t side, cham_dir_t dir, int M, int N, CHAMELEON_Complex64_t *A, int LDA, int K1, int K2, int *IPIV );
int CHAMELEON_zlatms( int M, int N, cham_dist_t idist, unsigned long long int seed, cham_sym_t sym, double *D, int mode, double cond, double dmax, CHAMELEON_Complex64_t *A, int LDA );
int CHAMELEON_zlatro( cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB );
int CHAMELEON_zlauum(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zplghe( double bump, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, unsigned long long int seed );
int CHAMELEON_zplgsy( CHAMELEON_Complex64_t bump, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, unsigned long long int seed );
int CHAMELEON_zplgtr( CHAMELEON_Complex64_t bump, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA, unsigned long long int seed );
int CHAMELEON_zplrnt( int M, int N, CHAMELEON_Complex64_t *A, int LDA, unsigned long long int seed );
int CHAMELEON_zplrnk( int M, int N, int K, CHAMELEON_Complex64_t *C, int LDC, unsigned long long int seedA, unsigned long long int seedB );
int CHAMELEON_zpoinv(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zposv(cham_uplo_t uplo, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zpotrf(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zsytrf(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zpotri(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zpotrimm(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zpotrs(cham_uplo_t uplo, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zsymm(cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zsyrk(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zsyr2k(cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC);
int CHAMELEON_zsysv(cham_uplo_t uplo, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zsytrs(cham_uplo_t uplo, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztpgqrt( int M, int N, int K, int L, CHAMELEON_Complex64_t *V1, int LDV1, CHAM_desc_t *descT1, CHAMELEON_Complex64_t *V2, int LDV2, CHAM_desc_t *descT2, CHAMELEON_Complex64_t *Q1, int LDQ1, CHAMELEON_Complex64_t *Q2, int LDQ2 );
int CHAMELEON_ztpqrt( int M, int N, int L, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAM_desc_t *descT );
int CHAMELEON_ztradd(cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztrmm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int N, int NRHS, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztrsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int N, int NRHS, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztrsmpl(int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descL, int *IPIV, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztrsmrv(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int N, int NRHS, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_ztrtri(cham_uplo_t uplo, cham_diag_t diag, int N, CHAMELEON_Complex64_t *A, int LDA);
int CHAMELEON_zunglq(int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zungqr(int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunmlq(cham_side_t side, cham_trans_t trans, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunmqr(cham_side_t side, cham_trans_t trans, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAMELEON_Complex64_t *B, int LDB);

//int CHAMELEON_zgecfi(int m, int n, CHAMELEON_Complex64_t *A, CHAMELEON_enum fin, int imb, int inb, CHAMELEON_enum fout, int omb, int onb);
//int CHAMELEON_zgetmi(int m, int n, CHAMELEON_Complex64_t *A, CHAMELEON_enum fin, int mb,  int nb);

/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
int CHAMELEON_zgeadd_Tile(cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B);
//int CHAMELEON_zgebrd_Tile(CHAM_desc_t *A, double *D, double *E, CHAM_desc_t *T);
//int CHAMELEON_zgecon_Tile(cham_normtype_t norm, CHAM_desc_t *A, double anorm, double *rcond);
//int CHAMELEON_zpocon_Tile(cham_uplo_t uplo, CHAM_desc_t *A, double anorm, double *rcond);
int CHAMELEON_zgelqf_Tile(CHAM_desc_t *A, CHAM_desc_t *T);
int CHAMELEON_zgelqs_Tile(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
int CHAMELEON_zgels_Tile(cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
double CHAMELEON_zgenm2_Tile( double tol, CHAM_desc_t *A );
int CHAMELEON_zgemm_Tile(cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C);
int CHAMELEON_zgepdf_qdwh_Tile( CHAM_desc_t *A, CHAM_desc_t *H, gepdf_info_t *info );
int CHAMELEON_zgepdf_qr_Tile( int doqr, int optid, const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB, CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1, CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2 );
int CHAMELEON_zgeqrf_Tile(CHAM_desc_t *A, CHAM_desc_t *T);
int CHAMELEON_zgeqrs_Tile(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
int CHAMELEON_zgesv_Tile(CHAM_desc_t *A, CHAM_ipiv_t *IPIV, CHAM_desc_t *B);
int CHAMELEON_zgesv_incpiv_Tile(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B);
int CHAMELEON_zgesv_nopiv_Tile(CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zgesvd_Tile(cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *A, double *S, CHAM_desc_t *T, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT);
//int CHAMELEON_zgetrf_Tile(CHAM_desc_t *A, int *IPIV);
int CHAMELEON_zgetrf_incpiv_Tile(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV);
int CHAMELEON_zgetrf_nopiv_Tile(CHAM_desc_t *A);
int CHAMELEON_zgetrf_Tile( CHAM_desc_t *A, CHAM_ipiv_t *IPIV );
//int CHAMELEON_zgetri_Tile(CHAM_desc_t *A, int *IPIV);
int CHAMELEON_zgetrs_Tile(cham_trans_t trans, CHAM_desc_t *A, CHAM_ipiv_t *IPIV, CHAM_desc_t *B);
int CHAMELEON_zgetrs_incpiv_Tile(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B);
int CHAMELEON_zgetrs_nopiv_Tile(CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zhemm_Tile(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C);
int CHAMELEON_zherk_Tile(cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C);
int CHAMELEON_zher2k_Tile(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta, CHAM_desc_t *C);
//int CHAMELEON_zheev_Tile(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *W, CHAM_desc_t *T, CHAMELEON_Complex64_t *Q, int LDQ);
int CHAMELEON_zheevd_Tile(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *W, CHAM_desc_t *T);
//int CHAMELEON_zhegv_Tile( CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, double *W, CHAM_desc_t *T, CHAM_desc_t *Q);
//int CHAMELEON_zhegvd_Tile(CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, double *W, CHAM_desc_t *T, CHAM_desc_t *Q);
//int CHAMELEON_zhegst_Tile(CHAMELEON_enum itype, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zhetrd_Tile(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *D, double *E, CHAM_desc_t *T, CHAMELEON_Complex64_t *Q, int LDQ);
int CHAMELEON_zlacpy_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
double CHAMELEON_zlange_Tile(cham_normtype_t norm, CHAM_desc_t *A);
double CHAMELEON_zlanhe_Tile(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A);
double CHAMELEON_zlansy_Tile(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A);
double CHAMELEON_zlantr_Tile(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A);
int CHAMELEON_zlascal_Tile(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A);
int CHAMELEON_zlaset_Tile(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *A);
int CHAMELEON_zlaswp_Tile( cham_side_t side, cham_dir_t dir, CHAM_desc_t *A, int K1, int K2, CHAM_ipiv_t *IPIV );
int CHAMELEON_zlatms_Tile( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym, double *D, int mode, double cond, double dmax, CHAM_desc_t *A );
int CHAMELEON_zlatro_Tile( cham_uplo_t uplo, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B );
int CHAMELEON_zlauum_Tile(cham_uplo_t uplo, CHAM_desc_t *A);
int CHAMELEON_zplghe_Tile(double bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed );
int CHAMELEON_zplgsy_Tile(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed );
int CHAMELEON_zplgtr_Tile(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed );
int CHAMELEON_zplrnt_Tile(CHAM_desc_t *A, unsigned long long int seed );
int CHAMELEON_zplrnk_Tile(int K, CHAM_desc_t *C, unsigned long long int seedA, unsigned long long int seedB );
int CHAMELEON_zpoinv_Tile(cham_uplo_t uplo, CHAM_desc_t *A);
int CHAMELEON_zposv_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zpotrf_Tile(cham_uplo_t uplo, CHAM_desc_t *A);
int CHAMELEON_zgered_Tile( cham_uplo_t uplo, double prec, CHAM_desc_t *A );
int CHAMELEON_zhered_Tile( cham_uplo_t uplo, double prec, CHAM_desc_t *A );
int CHAMELEON_zgerst_Tile( cham_uplo_t uplo, CHAM_desc_t *A );
int CHAMELEON_zsytrf_Tile(cham_uplo_t uplo, CHAM_desc_t *A);
int CHAMELEON_zpotri_Tile(cham_uplo_t uplo, CHAM_desc_t *A);
int CHAMELEON_zpotrimm_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C);
int CHAMELEON_zpotrs_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zsymm_Tile(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C);
int CHAMELEON_zsyrk_Tile(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *C);
int CHAMELEON_zsyr2k_Tile(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C);
int CHAMELEON_zsysv_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_zsytrs_Tile(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_ztpgqrt_Tile( int L, CHAM_desc_t *V1, CHAM_desc_t *T1, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1, CHAM_desc_t *Q2 );
int CHAMELEON_ztpqrt_Tile( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T );
int CHAMELEON_ztradd_Tile(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B);
int CHAMELEON_ztrmm_Tile(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_ztrsm_Tile(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_ztrsmpl_Tile(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B);
int CHAMELEON_ztrsmrv_Tile(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B);
int CHAMELEON_ztrtri_Tile(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A);
int CHAMELEON_zunglq_Tile(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
int CHAMELEON_zungqr_Tile(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
int CHAMELEON_zunmlq_Tile(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);
int CHAMELEON_zunmqr_Tile(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B);

/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
int CHAMELEON_zgeadd_Tile_Async(cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zgebrd_Tile_Async(CHAM_desc_t *A, double *D, double *E, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zgecon_Tile_Async(cham_normtype_t norm, CHAM_desc_t *A, double anorm, double *rcond, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zpocon_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, double anorm, double *rcond, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgelqf_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgelqs_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgels_Tile_Async(cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgenm2_Tile_Async( double tol, CHAM_desc_t *A, double *value, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zgemm_Tile_Async(cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgepdf_qdwh_Tile_Async( CHAM_desc_t *A, CHAM_desc_t *H, gepdf_info_t *info, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zgeqrf_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgeqrs_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgesv_Tile_Async(CHAM_desc_t *A, CHAM_ipiv_t *IPIV, CHAM_desc_t *B, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgesv_incpiv_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgesv_nopiv_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *B, void * ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgesvd_Tile_Async(cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *A, double *S, CHAM_desc_t *T, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *VT, int LDVT, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zgetrf_Tile_Async(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrf_incpiv_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrf_nopiv_Tile_Async(CHAM_desc_t *A, void * ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrf_Tile_Async( CHAM_desc_t *A, CHAM_ipiv_t *IPIV, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
//int CHAMELEON_zgetri_Tile_Async(CHAM_desc_t *A, int *IPIV, CHAM_desc_t *W, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrs_Tile_Async(cham_trans_t trans, CHAM_desc_t *A, CHAM_ipiv_t *IPIV, CHAM_desc_t *B, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrs_incpiv_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgetrs_nopiv_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zhemm_Tile_Async(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zherk_Tile_Async(cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zher2k_Tile_Async(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zheev_Tile_Async(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *W, CHAM_desc_t *T, CHAMELEON_Complex64_t *Q, int LDQ, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zheevd_Tile_Async(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *W, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zhegv_Tile_Async( CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, double *W, CHAM_desc_t *T, CHAM_desc_t *Q, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zhegvd_Tile_Async(CHAMELEON_enum itype, cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, double *W, CHAM_desc_t *T, CHAM_desc_t *Q, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zhegst_Tile_Async(CHAMELEON_enum itype, cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zhetrd_Tile_Async(cham_job_t jobz, cham_uplo_t uplo, CHAM_desc_t *A, double *D, double *E, CHAM_desc_t *T, CHAMELEON_Complex64_t *Q, int LDQ, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlacpy_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlange_Tile_Async(cham_normtype_t norm, CHAM_desc_t *A, double *value, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlanhe_Tile_Async(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *value, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlansy_Tile_Async(cham_normtype_t norm, cham_uplo_t uplo, CHAM_desc_t *A, double *value, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlantr_Tile_Async(cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, double *value, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlascal_Tile_Async(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlaset_Tile_Async(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zlaswp_Tile_Async( cham_side_t side, cham_dir_t dir, CHAM_desc_t *A, int K1, int K2, CHAM_ipiv_t *IPIV, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zlatms_Tile_Async( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym, double *D, int mode, double cond, double dmax, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zlatro_Tile_Async( cham_uplo_t uplo, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zlauum_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zplghe_Tile_Async(double bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zplgsy_Tile_Async(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zplgtr_Tile_Async(CHAMELEON_Complex64_t bump, cham_uplo_t uplo, CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zplrnt_Tile_Async(CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zplrnk_Tile_Async(int K, CHAM_desc_t *C, unsigned long long int seedA, unsigned long long int seedB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zpoinv_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zposv_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zpotrf_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgered_Tile_Async(cham_uplo_t uplo, double prec, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zhered_Tile_Async(cham_uplo_t uplo, double prec, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgerst_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_zsytrf_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zpotri_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zpotrimm_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zpotrs_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zsysv_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zsytrs_Tile_Async(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zsymm_Tile_Async(cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, void *ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zsyrk_Tile_Async(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zsyr2k_Tile_Async(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztpgqrt_Tile_Async( int L, CHAM_desc_t *V1, CHAM_desc_t *T1, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1, CHAM_desc_t *Q2, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_ztpqrt_Tile_Async( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
int CHAMELEON_ztradd_Tile_Async(cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAMELEON_Complex64_t beta, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztrmm_Tile_Async(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztrsm_Tile_Async(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztrsmpl_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *L, int *IPIV, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztrsmrv_Tile_Async(cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_ztrtri_Tile_Async(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunglq_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zungqr_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunmlq_Tile_Async(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunmqr_Tile_Async(cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

//int CHAMELEON_zgecfi_Async(int m, int n, CHAMELEON_Complex64_t *A, CHAMELEON_enum f_in, int imb, int inb, CHAMELEON_enum f_out, int omb, int onb, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
//int CHAMELEON_zgetmi_Async(int m, int n, CHAMELEON_Complex64_t *A, CHAMELEON_enum f_in, int mb,  int inb, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

/**
 *  Declarations of libhqr dependent functions.
 */
/**
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 */
int CHAMELEON_zgels_param(const libhqr_tree_t *qrtree, cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *descTT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgelqf_param(const libhqr_tree_t *qrtree, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *descTT);
int CHAMELEON_zgelqs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *descTT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zgeqrf_param(const libhqr_tree_t *qrtree, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *descTT);
int CHAMELEON_zgeqrs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *TT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunglq_param(const libhqr_tree_t *qrtree, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *TT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zungqr_param(const libhqr_tree_t *qrtree, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *descTT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunmlq_param(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAM_desc_t *TT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunmlq_param(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descT, CHAM_desc_t *TT, CHAMELEON_Complex64_t *B, int LDB);
int CHAMELEON_zunmqr_param(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descTS, CHAM_desc_t *TT, CHAMELEON_Complex64_t *B, int LDB);
/**
 *  Declarations of math functions (tile layout) - alphabetical order
 */
int CHAMELEON_zgels_param_Tile(const libhqr_tree_t *qrtree, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zgelqf_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT);
int CHAMELEON_zgelqs_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zgeqrf_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT);
int CHAMELEON_zgeqrs_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zunglq_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zungqr_param_Tile(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zunmlq_param_Tile(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
int CHAMELEON_zunmqr_param_Tile(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B);
/**
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 */
int CHAMELEON_zgels_param_Tile_Async(const libhqr_tree_t *qrtree, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgelqf_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgelqs_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgeqrf_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zgeqrs_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunglq_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zungqr_param_Tile_Async(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunmlq_param_Tile_Async(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);
int CHAMELEON_zunmqr_param_Tile_Async(const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

/**
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 */
void *CHAMELEON_zgemm_WS_Alloc( cham_trans_t transA, cham_trans_t transB, const CHAM_desc_t *A, const CHAM_desc_t *B, const CHAM_desc_t *C );
void  CHAMELEON_zgemm_WS_Free( void *ws );
void *CHAMELEON_zhemm_WS_Alloc( cham_side_t side, cham_uplo_t uplo, const CHAM_desc_t *A, const CHAM_desc_t *B, const CHAM_desc_t *C );
void  CHAMELEON_zhemm_WS_Free( void *ws );
void *CHAMELEON_zsymm_WS_Alloc( cham_side_t side, cham_uplo_t uplo, const CHAM_desc_t *A, const CHAM_desc_t *B, const CHAM_desc_t *C );
void  CHAMELEON_zsymm_WS_Free( void *ws );
void *CHAMELEON_zcesca_WS_Alloc( const CHAM_desc_t *A );
void  CHAMELEON_zcesca_WS_Free( void *ws );
void *CHAMELEON_zgram_WS_Alloc( const CHAM_desc_t *A );
void  CHAMELEON_zgram_WS_Free( void *ws );
void *CHAMELEON_zgetrf_WS_Alloc( const CHAM_desc_t *A );
void  CHAMELEON_zgetrf_WS_Free( void *ws );
void *CHAMELEON_zgetrf_nopiv_WS_Alloc( const CHAM_desc_t *A );
void  CHAMELEON_zgetrf_nopiv_WS_Free( void *ws );
void *CHAMELEON_zlaswp_WS_Alloc( cham_side_t side, const CHAM_desc_t *A );
void  CHAMELEON_zlaswp_WS_Free( void *ws );

int CHAMELEON_Alloc_Workspace_zgesv_incpiv(        int N, CHAM_desc_t **descL, int **IPIV, int p, int q);
int CHAMELEON_Alloc_Workspace_zgetrf_incpiv(int M, int N, CHAM_desc_t **descL, int **IPIV, int p, int q);

int CHAMELEON_Alloc_Workspace_zgebrd(int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgeev( int N,        CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgehrd(int N,        CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgelqf(int M, int N, CHAM_desc_t **T, int p, int q);
int CHAMELEON_Alloc_Workspace_zgels( int M, int N, CHAM_desc_t **T, int p, int q);
int CHAMELEON_Alloc_Workspace_zgeqrf(int M, int N, CHAM_desc_t **T, int p, int q);
int CHAMELEON_Alloc_Workspace_zgesvd(int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zheev( int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zheevd( int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zhegv( int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zhegvd( int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zhetrd(int M, int N, CHAM_desc_t **descT, int p, int q);

/**
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 */
int CHAMELEON_Alloc_Workspace_zgetri_Tile_Async(CHAM_desc_t *A, CHAM_desc_t *W);

/* Warning: Those functions are deprecated */
int CHAMELEON_Alloc_Workspace_zgelqf_Tile(int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgels_Tile( int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgeqrf_Tile(int M, int N, CHAM_desc_t **descT, int p, int q);
int CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile (int N, CHAM_desc_t **descL, int **IPIV, int p, int q);
int CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile(int N, CHAM_desc_t **descL, int **IPIV, int p, int q);

/**
 *  Auxiliary function prototypes
 */
int CHAMELEON_zLapack_to_Tile( CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descA ) __attribute__((deprecated("Please refer to CHAMELEON_zLap2Desc() instead")));
int CHAMELEON_zTile_to_Lapack( CHAM_desc_t *descA, CHAMELEON_Complex64_t *A, int LDA ) __attribute__((deprecated("Please refer to CHAMELEON_zDesc2Lap() instead")));
int CHAMELEON_zLap2Desc( cham_uplo_t uplo, CHAMELEON_Complex64_t *A, int LDA, CHAM_desc_t *descA );
int CHAMELEON_zDesc2Lap( cham_uplo_t uplo, CHAM_desc_t *descA, CHAMELEON_Complex64_t *A, int LDA );

/**
 * Centered-Scaled function prototypes
 */
int CHAMELEON_zcesca(int center, int scale, cham_store_t axis, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *SR, CHAMELEON_Complex64_t *SC );
int CHAMELEON_zcesca_Tile( int center, int scale, cham_store_t axis, CHAM_desc_t *A, CHAMELEON_Complex64_t *SR, CHAMELEON_Complex64_t *SC );
int CHAMELEON_zcesca_Tile_Async( int center, int scale, cham_store_t axis, CHAM_desc_t *A, void *user_ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );
/**
 * Gram function prototypes
 */
int CHAMELEON_zgram( cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA );
int CHAMELEON_zgram_Tile( cham_uplo_t uplo, CHAM_desc_t *A );
int CHAMELEON_zgram_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *A, void *user_ws, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

int CHAMELEON_zprint( FILE *file, const char *header, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA );
int CHAMELEON_zprint_Tile( FILE *file, const char *header, cham_uplo_t uplo, CHAM_desc_t *descA );

/**
 * Batch function prototypes - Tile interface
 */
int CHAMELEON_zgemm_batch_Tile( cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B, CHAMELEON_Complex64_t beta, CHAM_desc_t *C );
int CHAMELEON_zherk_batch_Tile( cham_uplo_t uplo, cham_trans_t trans, double alpha, CHAM_desc_t *A, double beta,  CHAM_desc_t *C );
int CHAMELEON_zplghe_batch_Tile( double bump, CHAM_desc_t *A, unsigned long long int seed );
int CHAMELEON_zpotrf_batch_Tile( cham_uplo_t uplo, CHAM_desc_t *A );
int CHAMELEON_ztrsm_batch_Tile( cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B );

END_C_DECLS

#endif /* _chameleon_z_h_ */
