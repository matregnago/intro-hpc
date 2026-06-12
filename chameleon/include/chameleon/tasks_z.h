/**
 *
 * @file chameleon_tasks_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t elementary tasks header
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
 * @author Romain Peressoni
 * @author Matthieu Kuhn
 * @author Ana Hourcau
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_tasks_z_h_
#define _chameleon_tasks_z_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB );
void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY );
void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                              cham_trans_t transA, cham_trans_t transB,
                              int m, int n, int k, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_zgered( const RUNTIME_option_t *options,
                         double threshold, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *Wnorm, int Wnm, int Wnn );
void INSERT_TASK_zgerst( const RUNTIME_option_t *options,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo );
void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zhemm_Astat( const RUNTIME_option_t *options,
                              cham_side_t side, cham_uplo_t uplo,
                              int m, int n, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zherfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb );
void INSERT_TASK_zlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn );
void INSERT_TASK_zlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn );
void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn );
void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t              dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *IPIV, int IPIVk,
                             const CHAM_desc_t *WA, int WAm, int WAn,
                             const CHAM_desc_t *A,  int Am,  int An );
void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 cham_side_t             side,
                                 cham_dir_t              dir,
                                 int                     m0,
                                 int                     m,
                                 int                     n,
                                 int                     k,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv, int ipivk,
                                 const CHAM_desc_t      *A,   int Am,   int An,
                                 const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                 const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                 void                  **clargs_ptr );
void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                       const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                       void                  **clargs_ptr );
void INSERT_TASK_zlaswp_get_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *A,    int Am,   int An,
                                     const CHAM_desc_t      *WAP,  int WAPm, int WAPn,
                                     void                  **clargs_ptr );
void INSERT_TASK_zlaswp_get_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WAP,  int WAPm, int WAPn,
                                           void                  **clargs_ptr );
void INSERT_TASK_zlaswp_set_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *A,    int Am,  int An,
                                     const CHAM_desc_t      *WA,   int WAm, int WAn,
                                     void                  **clargs_ptr );
void INSERT_TASK_zlaswp_set_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WA,   int WAm, int WAn,
                                           void                  **clargs_ptr );
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn );
void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo );
void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV );
void INSERT_TASK_zsymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsymm_Astat( const RUNTIME_option_t *options,
                              cham_side_t side, cham_uplo_t uplo,
                              int m, int n, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo );
void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,

                         int iinfo );
void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );

/**
 * Keep these insert_task for retro-compatibility
 */
static inline void
INSERT_TASK_ztslqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztplqt( options, m, n, 0, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_ztsqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztpqrt( options, m, n, 0, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_zttlqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztplqt( options, m, n, n, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_zttqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztpqrt( options, m, n, m, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_ztsmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, 0, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_ztsmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, 0, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_zttmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, n2, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_zttmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, m2, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

/**
 * Centered-Scaled and Gram prototype
 */
void INSERT_TASK_zgesum( const RUNTIME_option_t *options,
                        cham_store_t storev, int m, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SUMS, int SUMSm, int SUMSn );
void INSERT_TASK_zcesca( const RUNTIME_option_t *options,
                         int center, int scale, cham_store_t axis,
                         int m, int n, int mt, int nt,
                         const CHAM_desc_t *Gi, int Gim, int Gin,
                         const CHAM_desc_t *Gj, int Gjm, int Gjn,
                         const CHAM_desc_t *G, int Gm, int Gn,
                         const CHAM_desc_t *Di, int Dim, int Din,
                         const CHAM_desc_t *Dj, int Djm, int Djn,
                         CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An);

void RUNTIME_zgersum_set_methods( const CHAM_desc_t *A, int Am, int An );
void RUNTIME_zgersum_submit_tree( const RUNTIME_option_t *options,
                                  const CHAM_desc_t *A, int Am, int An );

/*
 * Tasks for LU factorization with partial pivoting
 */
void INSERT_TASK_zgetrf_nopiv_percol_diag( const RUNTIME_option_t *options,
                                           int m, int n, int k,
                                           const CHAM_desc_t *A, int Am, int An,
                                           const CHAM_desc_t *U, int Um, int Un,
                                           int iinfo );

void INSERT_TASK_zgetrf_nopiv_percol_trsm( const RUNTIME_option_t *options,
                                           int m, int n, int k,
                                           const CHAM_desc_t *A, int Am, int An,
                                           const CHAM_desc_t *U, int Um, int Un );

void INSERT_TASK_zgetrf_percol_diag( const RUNTIME_option_t *options,
                                     int m, int n, int h, int m0,
                                     CHAM_desc_t *A, int Am, int An,
                                     CHAM_ipiv_t       *ipiv,
                                     CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_percol_offdiag( const RUNTIME_option_t *options,
                                        int m, int n, int h, int m0,
                                        CHAM_desc_t *A, int Am, int An,
                                        CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib, int readUp,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv,
                                      CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib, int readUp,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_panel_offdiag_batched( const RUNTIME_option_t *options,
                                               int m, int n, int h, int m0,
                                               void *ws,
                                               CHAM_desc_t *A, int Am, int An,
                                               void **clargs_ptr,
                                               CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                     CHAM_desc_t *A, int An,
                                                     void **clargs_ptr,
                                                     CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_panel_blocked_batched( const RUNTIME_option_t *options,
                                               int m, int n, int h, int m0, int readUp,
                                               void *ws,
                                               CHAM_desc_t *A, int Am, int An,
                                               CHAM_desc_t *U, int Um, int Un,
                                               void **clargs_ptr,
                                               CHAM_ipiv_t *ipiv,
                                               CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_panel_blocked_batched_flush( const RUNTIME_option_t *options,
                                                     CHAM_desc_t *A, int An, int readUp,
                                                     CHAM_desc_t *U, int Um, int Un,
                                                     void **clargs_ptr,
                                                     CHAM_ipiv_t *ipiv,
                                                     CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int m, int n, int h, int ib,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_desc_pivot_t *pivot );

void INSERT_TASK_zgetrf_cpy_pivrow_in_Up( const RUNTIME_option_t *options,
                                          CHAM_desc_t            *Up,
                                          int                     Upm,
                                          int                     k,
                                          int                     h,
                                          int                     ib,
                                          int                     n,
                                          CHAM_desc_pivot_t      *pivot );

void INSERT_TASK_zlaswp_gemm( const RUNTIME_option_t *options,
                              cham_side_t             side,
                              cham_dir_t              dir,
                              int                     m0,
                              int                     m,
                              int                     n,
                              int                     k,
                              void                   *ws,
                              const CHAM_ipiv_t      *ipiv, int ipivk,
                              const CHAM_desc_t      *WA,   int WAm, int WAn,
                              const CHAM_desc_t      *A,    int Am,  int An,
                              const CHAM_desc_t      *B,    int Bm,  int Bn,
                              const CHAM_desc_t      *C,    int Cm,  int Cn,
                              void                  **clargs_ptr );

void INSERT_TASK_zlaswp_gemm_flush( const RUNTIME_option_t *options,
                                    cham_dir_t              dir,
                                    const CHAM_ipiv_t      *ipiv, int ipivk,
                                    const CHAM_desc_t      *WA,   int WAm, int WAn,
                                    const CHAM_desc_t      *B,    int Bm,  int Bn,
                                    void                  **clargs_ptr );

#if defined(CHAMELEON_USE_MPI)
/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Perfoms an allreduce operation on the pivot k during the panel
 *  factorization of the LU factorization with partial pivoting.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] A
 *          The descriptor of the matrix A.

 * @param[in] pivot
 *          The pivot structure that contains the informations for the LU
 *          factorization with partial pivoting.
*
 * @param[in] k
 *          The iteration k of the panel factorizatio).
 *
 * @param[in] h
 *          The iteration h of the column of the panel factorization.
 *
 * @param[in] n
 *          The number of columns of the row.
 *
 * @param[in] ws
 *          The workspace to handle the data in the LU factorization with
 *          partial pivoting.
 *
 *******************************************************************************
 */
void INSERT_TASK_zpivot_allreduce( const RUNTIME_option_t *options,
                                   CHAM_desc_t            *A,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     k,
                                   int                     h,
                                   int                     n,
                                   void                   *ws );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Perfoms an allreduce operation on the tile
 *  U(Um, Un) according to the permutation ipiv. This task is used in the LU
 *  factorization with partial pivoting.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[inout] U
 *          The descriptor of the worskpace used for the permutation in the LU
 *          factorization with partial pivoting.
 *
 * @param[in] Um
 *          The row index of the tile used in U.
 *
 * @param[in] Un
 *          The column index of the tile used in U.
 *
 * @param[in] ipiv
 *          The pivot structure that contains the informations for the LU
 *          factorization with partial pivoting.
 *
 * @param[in] ipivk
 *          The index of the permutation.
 *
 * @param[in] k
 *          The number of rows in the tile U(Um, Un).
 *
 * @param[in] n
 *          The number of columns in the tile U(Um, Un).
 *
 * @param[in] ws
 *          The workspace to handle the data in the LU factorization with
 *          partial pivoting.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
                                  cham_dir_t              dir,
                                  const CHAM_desc_t      *A,
                                  CHAM_desc_t            *U,
                                  int                     Um,
                                  int                     Un,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     ipivk,
                                  int                     k,
                                  int                     n,
                                  void                   *ws );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Perfoms an reduce operation on the tile
 *  U(Um, Un) according to the permutation ipiv. This task is used in the
 *  LASWP.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] Am
 *          The index of the row of A to reduce.
 *
 * @param[in] An
 *          The index of the column of A to reduce.
 *
 * @param[in] ipiv
 *          The descriptor of array of pivot ipiv.
 *
 * @param[in] ipivk
 *          The index of the current tile of the pivot ipiv.
 *
 * @param[in] Wu
 *          The workspace to handle the data in the LU factorization with
 *          partial pivoting.
 *
 * @param[in] Wum
 *          The row index of Wu.
 *
 * @param[in] Wun
 *          The column index of Wu.
 *
 * @param[in] ws
 *          The workspace to handle the data in the LU factorization with
 *          partial pivoting.
 *
 * @param[in] Wm
 *          The row index of Ws.
 *
 * @param[in] Wn
 *          The column index of Ws.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_reduce( const RUNTIME_option_t *options,
                               cham_dir_t              dir,
                               const CHAM_desc_t      *A,
                               int                     Am,
                               int                     An,
                               CHAM_ipiv_t            *ipiv,
                               int                     ipivk,
                               const CHAM_desc_t      *Wu,
                               int                     Wum,
                               int                     Wun,
                               void                   *ws,
                               int                     Wm,
                               int                     Wn );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Sends the tile A(Am, An) to the processus
 *  involved in the permutation. This task is used in the LU factorization with
 *  partial pivoting.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] Am
 *          The row index of the tile used in A.
 *
 * @param[in] An
 *          The column index of the tile used in A.
 *
 * @param[in] myrank
 *          The rank of the current process.
 *
 * @param[in] np
 *          The number of processus involved in the permutation.
 *
 * @param[in] proc_involved
 *          The list of the processus involved in the permutation.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                         CHAM_desc_t            *A,
                                         int                     Am,
                                         int                     An,
                                         int                     myrank,
                                         int                     np,
                                         int                    *proc_involved );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief - Sends the permutation ipivk to the
 *  processus involved in the permutation. This task is used in the LU
 *  factorization with partial pivoting.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order
 *          = ChamDirBackward: Reverse order
 *
 * @param[in] ipiv
 *          The pivot structure that contains the informations for the LU
 *          factorization with partial pivoting.
 *
 * @param[in] ipivk
 *          The index of the permutation.
 *
 * @param[in] myrank
 *          The rank of the current process.
 *
 * @param[in] np
 *          The number of processus involved in the permutation.
 *
 * @param[in] proc_involved
 *          The list of the processus involved in the permutation.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                            cham_dir_t              dir,
                                            CHAM_ipiv_t            *ipiv,
                                            int                     ipivk,
                                            int                     myrank,
                                            int                     np,
                                            int                    *proc_involved );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Sends the inverse permutation ipivk
 *  to the processus involved in the permutation. This task is used in the LU
 *  factorization with partial pivoting.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order
 *          = ChamDirBackward: Reverse order
 *
 * @param[in] ipiv
 *          The pivot structure that contains the informations for the LU
 *          factorization with partial pivoting.
 *
 * @param[in] ipivk
 *          The index of the permutation.
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] k
 *          The index of the panel factorized.
 *
 * @param[in] n
 *          The index of the panel to permute.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_allreduce_send_invp_row( const RUNTIME_option_t *options,
                                                cham_dir_t              dir,
                                                CHAM_ipiv_t            *ipiv,
                                                int                     ipivk,
                                                const CHAM_desc_t      *A,
                                                int                     k,
                                                int                     n );

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Sends the inverse permutation ipivk
 *  to the processus involved in the permutation. This task is used in LASWP with
 *  columns permutation.
 *
 *******************************************************************************
 *
 * @param[in] options
 *          The runtime options data structure to pass through all insert_task calls.
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order
 *          = ChamDirBackward: Reverse order
 *
 * @param[in] ipiv
 *          The pivot structure that contains the informations for the LU
 *          factorization with partial pivoting.
 *
 * @param[in] ipivk
 *          The index of the permutation.
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] m
 *          The index of the panel to permute.
 *
 * @param[in] k
 *          The index of the panel factorized.
 *
 *******************************************************************************
 */
void INSERT_TASK_zperm_allreduce_send_invp_col( const RUNTIME_option_t *options,
                                                cham_dir_t              dir,
                                                CHAM_ipiv_t            *ipiv,
                                                int                     ipivk,
                                                const CHAM_desc_t      *A,
                                                int                     m,
                                                int                     k );

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws,    int Wm,     int Wn,
                             const CHAM_desc_t *tileA, int tileAm, int tileAn );

#endif /* defined(CHAMELEON_USE_MPI) */

#endif /* _chameleon_tasks_z_h_ */

