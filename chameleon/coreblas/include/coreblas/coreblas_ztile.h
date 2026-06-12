/**
 *
 * @file coreblas_ztile.h
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800 ), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CPU kernel CHAM_tile_t interface
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _coreblas_ztile_h_
#define _coreblas_ztile_h_

#if defined(PRECISION_z) || defined(PRECISION_c)
void TCORE_dlag2z( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B );
#endif
void TCORE_dzasum( cham_store_t storev, cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, double *work );
int  TCORE_zaxpy( int M, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, int incA, CHAM_tile_t *B, int incB );
int  TCORE_zgeadd( cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, CHAMELEON_Complex64_t beta, CHAM_tile_t *B );
int  TCORE_zgelqt( int M, int N, int IB, CHAM_tile_t *A, CHAM_tile_t *T, CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK );
void TCORE_zgemv( cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *x, int incx, CHAMELEON_Complex64_t beta, CHAM_tile_t *y, int incy );
void TCORE_zgemm( cham_trans_t transA, cham_trans_t transB, int M, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, CHAMELEON_Complex64_t beta, CHAM_tile_t *C );
int  TCORE_zgeqrt( int M, int N, int IB, CHAM_tile_t *A, CHAM_tile_t *T, CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK );
int  TCORE_zgessm( int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t *L, CHAM_tile_t *A );
int  TCORE_zgessq( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq );
int  TCORE_zgetrf( int M, int N, CHAM_tile_t *A, int *IPIV, int *INFO );
int  TCORE_zgetrf_incpiv( int M, int N, int IB, CHAM_tile_t *A, int *IPIV, int *INFO );
int  TCORE_zgetrf_nopiv( int M, int N, int IB, CHAM_tile_t *A, int *INFO );
void TCORE_zhe2ge( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B );
#if defined(PRECISION_z) || defined(PRECISION_c)
void TCORE_zhemm( cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, CHAMELEON_Complex64_t beta, CHAM_tile_t *C );
void TCORE_zherk( cham_uplo_t uplo, cham_trans_t trans, int N, int K, double alpha, const CHAM_tile_t *A, double beta, CHAM_tile_t *C );
void TCORE_zher2k( cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, double beta, CHAM_tile_t *C );
#endif
int  TCORE_zherfb( cham_uplo_t uplo, int N, int K, int IB, int NB, const CHAM_tile_t *A, const CHAM_tile_t *T, CHAM_tile_t *C, CHAMELEON_Complex64_t *WORK, int ldwork );
#if defined(PRECISION_z) || defined(PRECISION_c)
int  TCORE_zhessq( cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq );
#endif
void TCORE_zlacpy( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B );
void TCORE_zlacpyx( cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t *A, int LDA, int displB, CHAM_tile_t *B, int LDB );
void TCORE_zlange( cham_normtype_t norm, int M, int N, const CHAM_tile_t *A, double *work, double *normA );
#if defined(PRECISION_z) || defined(PRECISION_c)
void TCORE_zlanhe( cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t *A, double *work, double *normA );
#endif
void TCORE_zlansy( cham_normtype_t norm, cham_uplo_t uplo, int N, const CHAM_tile_t *A, double *work, double *normA );
void TCORE_zlantr( cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t *A, double *work, double *normA );
int  TCORE_zlascal( cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A );
void TCORE_zlaset( cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta, CHAM_tile_t *A );
void TCORE_zlaset2( cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A );
int TCORE_zlaswp_get( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *perm );
int TCORE_zlaswp_get_idx( cham_side_t side, int m0, int m, int n, int k,
                          const CHAM_tile_t *A, CHAM_tile_t *B,
                          int perm_m, int perm_mt, const int *perm_idx );
int TCORE_zlaswp_set( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *invp );
int TCORE_zlaswp_set_idx( cham_side_t side, int m0, int m, int n, int k,
                          const CHAM_tile_t *A, CHAM_tile_t *B,
                          int invp_m, int invp_mt, const int *invp_idx );
int  TCORE_zlatro( cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B );
void TCORE_zlauum( cham_uplo_t uplo, int N, CHAM_tile_t *A );
#if defined(PRECISION_z) || defined(PRECISION_c)
void TCORE_zplghe( double bump, int m, int n, CHAM_tile_t *tileA, int bigM, int m0, int n0, unsigned long long int seed );
#endif
void TCORE_zplgsy( CHAMELEON_Complex64_t bump, int m, int n, CHAM_tile_t *tileA, int bigM, int m0, int n0, unsigned long long int seed );
void TCORE_zplrnt( int m, int n, CHAM_tile_t *tileA, int bigM, int m0, int n0, unsigned long long int seed );
void TCORE_zpotrf( cham_uplo_t uplo, int n, CHAM_tile_t *A, int *INFO );
int  TCORE_zssssm( int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *L1, const CHAM_tile_t *L2, const int *IPIV );
void TCORE_zsymm( cham_side_t side, cham_uplo_t uplo, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, CHAMELEON_Complex64_t beta, CHAM_tile_t *C );
void TCORE_zsyrk( cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, CHAMELEON_Complex64_t beta, CHAM_tile_t *C );
void TCORE_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, CHAMELEON_Complex64_t beta, CHAM_tile_t *C );
int  TCORE_zsyssq( cham_store_t storev, cham_uplo_t uplo, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq );
int  TCORE_zsytf2_nopiv( cham_uplo_t uplo, int n, CHAM_tile_t *A );
int  TCORE_ztplqt( int M, int N, int L, int IB, CHAM_tile_t *A, CHAM_tile_t *B, CHAM_tile_t *T, CHAMELEON_Complex64_t *WORK );
int  TCORE_ztpmlqt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *A, CHAM_tile_t *B, CHAMELEON_Complex64_t *WORK );
int  TCORE_ztpmqrt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *A, CHAM_tile_t *B, CHAMELEON_Complex64_t *WORK );
int  TCORE_ztpqrt( int M, int N, int L, int IB, CHAM_tile_t *A, CHAM_tile_t *B, CHAM_tile_t *T, CHAMELEON_Complex64_t *WORK );
int  TCORE_ztradd( cham_uplo_t uplo, cham_trans_t trans, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, CHAMELEON_Complex64_t beta, CHAM_tile_t *B );
void TCORE_ztrasm( cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t *A, double *work );
void TCORE_ztrmm( cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, CHAM_tile_t *B );
void TCORE_ztrsm( cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int M, int N, CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A, CHAM_tile_t *B );
int  TCORE_ztrssq( cham_uplo_t uplo, cham_diag_t diag, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq );
void TCORE_ztrtri( cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t *A, int *info );
int  TCORE_ztsmlq_hetra1( cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAMELEON_Complex64_t *WORK, int ldwork );
int  TCORE_ztsmqr_hetra1( cham_side_t side, cham_trans_t trans, int m1, int n1, int m2, int n2, int k, int ib, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAMELEON_Complex64_t *WORK, int ldwork );
int  TCORE_ztstrf( int M, int N, int IB, int NB, CHAM_tile_t *U, CHAM_tile_t *A, CHAM_tile_t *L, int *IPIV, CHAMELEON_Complex64_t *WORK, int LDWORK, int *INFO );
int  TCORE_zunmlq( cham_side_t side, cham_trans_t trans, int M, int N, int IB, int K, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *C, CHAMELEON_Complex64_t *WORK, int LDWORK );
int  TCORE_zunmqr( cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *C, CHAMELEON_Complex64_t *WORK, int LDWORK );
int  TCORE_zgesum( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sum );
int  TCORE_zcesca( int center, int scale, cham_store_t axis, int M, int N, int Mt, int Nt, const CHAM_tile_t *Gi, const CHAM_tile_t *Gj, const CHAM_tile_t *G, const CHAM_tile_t *Di, const CHAM_tile_t *Dj, CHAM_tile_t *A );
int  TCORE_zgram( cham_uplo_t uplo, int M, int N, int Mt, int Nt, const CHAM_tile_t *Di, const CHAM_tile_t *Dj, const CHAM_tile_t *D, CHAM_tile_t *A );
void TCORE_zprint( FILE *file, const char *header, cham_uplo_t uplo, int M, int N, int Am, int An, const CHAM_tile_t *A );

#endif /* _coreblas_ztile_h_ */
