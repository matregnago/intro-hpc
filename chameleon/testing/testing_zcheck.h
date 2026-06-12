/**
 *
 * @file testing_zcheck.h
 *
 * @copyright 2019-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings header
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */

#ifndef _testing_zcheck_h_
#define _testing_zcheck_h_

#include "testings.h"
#include <math.h>

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif

#define CHECK_TRMM 3
#define CHECK_TRSM 4

#if defined(CHAMELEON_SIMULATION)

/* Auxiliary check */
static inline int check_zmatrices_std ( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB ) { return 0; }
static inline int check_zmatrices     ( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descB ) { return 0; }
static inline int check_znorm_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_normtype_t norm_type, cham_uplo_t uplo,
                                        cham_diag_t diag, double norm_cham, int M, int N, CHAMELEON_Complex64_t *A, int LDA ) { return 0; }
static inline int check_znorm         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_normtype_t norm_type, cham_uplo_t uplo,
                                        cham_diag_t diag, double norm_cham, CHAM_desc_t *descA ) { return 0; }
static inline int check_zsum_std      ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A,
                                        int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB ) { return 0; }
static inline int check_zsum          ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                                        CHAMELEON_Complex64_t beta, CHAM_desc_t *descBref, CHAM_desc_t *descBcham ) { return 0; }
static inline int check_zscale_std    ( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA ) { return 0; }
static inline int check_zscale        ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descAinit, CHAM_desc_t *descA ) { return 0; }
static inline int check_zrankk_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA ) { return 0; }
static inline int check_zrankk        ( run_arg_list_t *args, int K, CHAM_desc_t *descA ) { return 0; }
static inline int check_ztranspose_std( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LBD ) { return 0; }
static inline int check_ztranspose    ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descB ) { return 0; }
static inline int check_zset_std      ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, int M, int N, CHAMELEON_Complex64_t *A, int LDA ) { return 0; }
static inline int check_zset          ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *descA ) { return 0; }

/* Blas check */
static inline int check_zgemm_std     ( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA,
                                        CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC )  { return 0; }
static inline int check_zgemm         ( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                                        CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC ) { return 0; }
static inline int check_zsymm_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                                        CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC ) { return 0; }
static inline int check_zsymm         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                                        CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC ) { return 0; }
static inline int check_zsyrk_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                                        CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC ) { return 0; }
static inline int check_zsyrk         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                                        CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC ) { return 0; }
static inline int check_ztrmm_std     ( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha,
                                        CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *B, int LDB ) { return 0; }
static inline int check_ztrmm         ( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descBref, CHAM_desc_t *descB ) { return 0; }

/* Factorisation (Cholesky) check */
static inline int check_zlauum        ( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descAAt ) { return 0; }
static inline int check_zlauum_std    ( run_arg_list_t *args, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AAt, int LDA ) { return 0; }
static inline int check_zxxtrf        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, CHAM_desc_t *descA1, CHAM_desc_t *descA2 ) { return 0; }
static inline int check_zxxtrf_std    ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *LU, int LDA ) { return 0; }
static inline int check_zsolve        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo,
                                        CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB ) { return 0; }
static inline int check_zsolve_std    ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo, int N, int NRHS,
                                        CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, CHAMELEON_Complex64_t *B, int LDB ) { return 0; }
static inline int check_ztrtri        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_diag_t diag,
                                        CHAM_desc_t *descA, CHAM_desc_t *descAi ) { return 0; }
static inline int check_ztrtri_std    ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                                        int N, CHAMELEON_Complex64_t *A0, CHAMELEON_Complex64_t *Ai, int LDA ) { return 0; }

/* QR/LQ factorisation check */
static inline int check_zortho        ( run_arg_list_t *args, CHAM_desc_t *descQ ) { return 0; }
static inline int check_zortho_std    ( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *Q, int LDQ ) { return 0; }
static inline int check_zgeqrf        ( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ ) { return 0; }
static inline int check_zgeqrf_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ ) { return 0; }
static inline int check_zgelqf        ( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ ) { return 0; }
static inline int check_zgelqf_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ ) { return 0; }
static inline int check_zgels         ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB ) { return 0; }
static inline int check_zgels_std     ( run_arg_list_t *args, cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, int LDX, CHAMELEON_Complex64_t *B, int LDB ) { return 0; }
static inline int check_zgeqrs        ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR ) { return 0; }
static inline int check_zgelqs        ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR ) { return 0; }
static inline int check_zqc           ( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, CHAM_desc_t *descC, CHAM_desc_t *descQ, CHAM_desc_t *descCC ) { return 0; }
static inline int check_zqc_std       ( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, int M, int N,
                                        CHAMELEON_Complex64_t *C, CHAMELEON_Complex64_t *CC, int LDC, CHAMELEON_Complex64_t *Q, int LDQ ) { return 0; }

/* SVD check */
static inline int check_zgesvd_std    ( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, int M, int N, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA,
                                        double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt ) { return 0; }
static inline int check_zgesvd        ( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *descAinit, CHAM_desc_t *descA,
                                        double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt ) { return 0; }

/* Polar decomposition check */
static inline int check_zgepdf_qr     ( run_arg_list_t *args, CHAM_desc_t *descA1, CHAM_desc_t *descA2,
                                        CHAM_desc_t *descQ1, CHAM_desc_t *descQ2, CHAM_desc_t *descAF1 ) { return 0; }
static inline int check_zxxpd         ( run_arg_list_t *args,
                                        CHAM_desc_t *descA, CHAM_desc_t *descU, CHAM_desc_t *descH ) { return 0; }
static inline int check_zxxpd_std     ( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *A,
                                        CHAMELEON_Complex64_t *U, int LDA, CHAMELEON_Complex64_t *H, int LDH ) { return 0; }

#else /* !defined(CHAMELEON_SIMULATION) */

/* Auxiliary check */
int check_zmatrices_std ( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB );
int check_zmatrices     ( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descB );
int check_znorm_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_normtype_t norm_type, cham_uplo_t uplo,
                          cham_diag_t diag, double norm_cham, int M, int N, CHAMELEON_Complex64_t *A, int LDA );
int check_znorm         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_normtype_t norm_type, cham_uplo_t uplo,
                          cham_diag_t diag, double norm_cham, CHAM_desc_t *descA );
int check_zsum_std      ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A,
                          int LDA, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *Bcham, int LDB );
int check_zsum          ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAMELEON_Complex64_t beta, CHAM_desc_t *descBref, CHAM_desc_t *descBcham );
int check_zscale_std    ( run_arg_list_t *args, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA );
int check_zscale        ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descAinit, CHAM_desc_t *descA );
int check_zrankk_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA );
int check_zrankk        ( run_arg_list_t *args, int K, CHAM_desc_t *descA );
int check_ztranspose_std( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LBD );
int check_ztranspose    ( run_arg_list_t *args, cham_uplo_t uplo, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descB );
int check_zset_std      ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, int M, int N, CHAMELEON_Complex64_t *A, int LDA );
int check_zset          ( run_arg_list_t *args, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_desc_t *descA );

/* Blas check */
int check_zgemm_std     ( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA,
                          CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC );
int check_zgemm         ( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_zsymm_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                          CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC );
int check_zsymm         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_side_t side, cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                          CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_zsyrk_std     ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                          CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC );
int check_zsyrk         ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_trans_t trans, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                          CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC );
int check_ztrmm_std     ( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha,
                          CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *B, int LDB );
int check_ztrmm         ( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descBref, CHAM_desc_t *descB );

/* Factorisation (Cholesky) check */
int check_zlauum        ( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descAAt );
int check_zlauum_std    ( run_arg_list_t *args, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AAt, int LDA );
int check_zxxtrf        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, CHAM_desc_t *descA1, CHAM_desc_t *descA2 );
int check_zxxtrf_std    ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *LU, int LDA );
int check_zsolve        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo,
                          CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB );
int check_zsolve_std    ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_trans_t trans, cham_uplo_t uplo, int N, int NRHS,
                          CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, CHAMELEON_Complex64_t *B, int LDB );
int check_ztrtri        ( run_arg_list_t *args, cham_mtxtype_t mtxtype, cham_uplo_t uplo, cham_diag_t diag,
                          CHAM_desc_t *descA, CHAM_desc_t *descAi );
int check_ztrtri_std    ( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                          int N, CHAMELEON_Complex64_t *A0, CHAMELEON_Complex64_t *Ai, int LDA );

/* QR/LQ factorisation check */
int check_zortho        ( run_arg_list_t *args, CHAM_desc_t *descQ );
int check_zortho_std    ( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *Q, int LDQ );
int check_zgeqrf        ( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ );
int check_zgeqrf_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ );
int check_zgelqf        ( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ );
int check_zgelqf_std    ( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ );
int check_zgels         ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB );
int check_zgels_std     ( run_arg_list_t *args, cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, int LDX, CHAMELEON_Complex64_t *B, int LDB );
int check_zgeqrs        ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR );
int check_zgelqs        ( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR );
int check_zqc           ( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, CHAM_desc_t *descC, CHAM_desc_t *descQ, CHAM_desc_t *descCC );
int check_zqc_std       ( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, int M, int N,
                          CHAMELEON_Complex64_t *C, CHAMELEON_Complex64_t *CC, int LDC, CHAMELEON_Complex64_t *Q, int LDQ );

/* SVD check */
int check_zgesvd_std    ( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, int M, int N, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA,
                          double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt );
int check_zgesvd        ( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *descAinit, CHAM_desc_t *descA,
                          double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt );
/* Polar decomposition check */
int check_zgepdf_qr     ( run_arg_list_t *args, CHAM_desc_t *descA1, CHAM_desc_t *descA2,
                          CHAM_desc_t *descQ1, CHAM_desc_t *descQ2, CHAM_desc_t *descAF1 );
int check_zxxpd         ( run_arg_list_t *args,
                          CHAM_desc_t *descA, CHAM_desc_t *descU, CHAM_desc_t *descH );
int check_zxxpd_std     ( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *A,
                          CHAMELEON_Complex64_t *U, int LDA, CHAMELEON_Complex64_t *H, int LDH );

#endif /* defined(CHAMELEON_SIMULATION) */

#endif /* _testing_zcheck_h_ */
