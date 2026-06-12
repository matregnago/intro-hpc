/**
 *
 * @file core_ztile_empty.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CPU kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Abel Calluaud
 * @author Matteo Marcos
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/coreblas_ztile.h"

#if defined( CHAMELEON_USE_HMATOSS )
#include "coreblas/hmat.h"
#endif

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_dlag2z( __attribute__((unused)) cham_uplo_t uplo,
              __attribute__((unused)) int M,
              __attribute__((unused)) int N,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    CORE_dlag2z( uplo, M, N, A->mat, A->ld, B->mat, B->ld );
}
#endif

void
TCORE_dzasum( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_dzasum( storev, uplo, M, N, CHAM_tile_get_ptr( A ), A->ld, work );
}

int
TCORE_zaxpy( int                   M,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             int                   incA,
             CHAM_tile_t *         B,
             int                   incB )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgeadd( __attribute__((unused)) cham_trans_t          trans,
              __attribute__((unused)) int                   M,
              __attribute__((unused)) int                   N,
              __attribute__((unused)) CHAMELEON_Complex64_t alpha,
              __attribute__((unused)) const CHAM_tile_t *   A,
              __attribute__((unused)) CHAMELEON_Complex64_t beta,
              __attribute__((unused)) CHAM_tile_t *         B )
{
    coreblas_kernel_trace( A, B );
    return 0;
}

int
TCORE_zgelqt( __attribute__((unused)) int                    M,
              __attribute__((unused)) int                    N,
              __attribute__((unused)) int                    IB,
              __attribute__((unused)) CHAM_tile_t *          A,
              __attribute__((unused)) CHAM_tile_t *          T,
              __attribute__((unused)) CHAMELEON_Complex64_t *TAU,
              __attribute__((unused)) CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( A, T );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_zgemv( cham_trans_t trans, int M, int N,
             CHAMELEON_Complex64_t alpha, const CHAM_tile_t *A,
                                          const CHAM_tile_t *x, int incX,
             CHAMELEON_Complex64_t beta,        CHAM_tile_t *y, int incY )
{
    coreblas_kernel_trace( A, x, y );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( x->format & CHAMELEON_TILE_FULLRANK );
    assert( y->format & CHAMELEON_TILE_FULLRANK );
}

void
TCORE_zgemm( cham_trans_t          transA,
             cham_trans_t          transB,
             int                   M,
             int                   N,
             int                   K,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    coreblas_kernel_trace( A, B, C );
}

int
TCORE_zgeqrt( int                    M,
              int                    N,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *TAU,
              CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( A, T );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgessm( int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t *L, CHAM_tile_t *A )
{
    coreblas_kernel_trace( L, A );
    assert( L->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgessq( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq )
{
    coreblas_kernel_trace( A, sclssq );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sclssq->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgetrf( int M, int N, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    coreblas_kernel_trace( A );
    return 0;
}

int
TCORE_zgetrf_incpiv( int M, int N, int IB, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    coreblas_kernel_trace( A );
    return 0;
}

int
TCORE_zgetrf_nopiv( int M, int N, int IB, CHAM_tile_t *A, int *INFO )
{
    coreblas_kernel_trace( A );
    return 0;
}

void
TCORE_zhe2ge( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zhemm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zherk( cham_uplo_t        uplo,
             cham_trans_t       trans,
             int                N,
             int                K,
             double             alpha,
             const CHAM_tile_t *A,
             double             beta,
             CHAM_tile_t *      C )
{
    coreblas_kernel_trace( A, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zher2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   N,
              int                   K,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              const CHAM_tile_t *   B,
              double                beta,
              CHAM_tile_t *         C )
{
    coreblas_kernel_trace( A, B, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}
#endif

int
TCORE_zherfb( cham_uplo_t            uplo,
              int                    N,
              int                    K,
              int                    IB,
              int                    NB,
              const CHAM_tile_t *    A,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    ldwork )
{
    coreblas_kernel_trace( A, T, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCORE_zhessq( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    coreblas_kernel_trace( A, sclssq );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sclssq->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}
#endif

void
TCORE_zlacpy( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    return;
}

void
TCORE_zlacpyx( cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t *A, int LDA, int displB, CHAM_tile_t *B, int LDB )
{
    return;
}

void
TCORE_zlange( __attribute__((unused)) cham_normtype_t    norm,
              __attribute__((unused)) int                M,
              __attribute__((unused)) int                N,
              const CHAM_tile_t *A,
              __attribute__((unused)) double *           work,
              __attribute__((unused)) double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zlanhe( __attribute__((unused)) cham_normtype_t    norm,
              __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) int                N,
              const CHAM_tile_t *A,
              __attribute__((unused)) double *           work,
              __attribute__((unused)) double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}
#endif

void
TCORE_zlansy( __attribute__((unused)) cham_normtype_t    norm,
              __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) int                N,
              const CHAM_tile_t *A,
              __attribute__((unused)) double *           work,
              __attribute__((unused)) double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zlantr( __attribute__((unused)) cham_normtype_t    norm,
              __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) cham_diag_t        diag,
              __attribute__((unused)) int                M,
              __attribute__((unused)) int                N,
              const CHAM_tile_t *A,
              __attribute__((unused)) double *           work,
              __attribute__((unused)) double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

int
TCORE_zlascal( cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_zlaset( cham_uplo_t           uplo,
              int                   n1,
              int                   n2,
              CHAMELEON_Complex64_t alpha,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zlaset2( cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

int
TCORE_zlaswp_get( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *perm )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zlaswp_get_idx( __attribute__((unused)) cham_side_t side,
                      __attribute__((unused)) int         m0,
                      __attribute__((unused)) int         m,
                      __attribute__((unused)) int         n,
                      __attribute__((unused)) int         k,
                      const CHAM_tile_t                  *A,
                      CHAM_tile_t                        *B,
                      __attribute__((unused)) int        perm_m,
                      __attribute__((unused)) int        perm_mt,
                      __attribute__((unused)) const int *perm_idx )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zlaswp_set( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *invp )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zlaswp_set_idx( __attribute__((unused)) cham_side_t side,
                      __attribute__((unused)) int         m0,
                      __attribute__((unused)) int         m,
                      __attribute__((unused)) int         n,
                      __attribute__((unused)) int         k,
                      const CHAM_tile_t                  *A,
                      CHAM_tile_t                        *B,
                      __attribute__((unused)) int        invp_m,
                      __attribute__((unused)) int        invp_mt,
                      __attribute__((unused)) const int *invp_idx )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zlatro( cham_uplo_t        uplo,
              cham_trans_t       trans,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_zlauum( cham_uplo_t uplo, int N, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zplghe( double                 bump,
              int                    m,
              int                    n,
              CHAM_tile_t *          A,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}
#endif

void
TCORE_zplgsy( CHAMELEON_Complex64_t  bump,
              int                    m,
              int                    n,
              CHAM_tile_t *          A,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zplrnt( int                    m,
              int                    n,
              CHAM_tile_t *          A,
              int                    bigM,
              int                    m0,
              int                    n0,
              unsigned long long int seed )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zpotrf( cham_uplo_t uplo, int n, CHAM_tile_t *A, int *INFO )
{
    coreblas_kernel_trace( A );
    return;
}

int
TCORE_zssssm( int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *L1,
              const CHAM_tile_t *L2,
              const int *        IPIV )
{
    coreblas_kernel_trace( A1, A2, L1, L2 );
    assert( A1->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A2->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( L1->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( L2->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_zsymm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             const CHAM_tile_t *   B,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    coreblas_kernel_trace( A, B, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zsyrk( cham_uplo_t           uplo,
             cham_trans_t          trans,
             int                   N,
             int                   K,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAMELEON_Complex64_t beta,
             CHAM_tile_t *         C )
{
    coreblas_kernel_trace( A, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_zsyr2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   N,
              int                   K,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              const CHAM_tile_t *   B,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         C )
{
    coreblas_kernel_trace( A, B, C );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

int
TCORE_zsyssq( cham_store_t       storev,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    coreblas_kernel_trace( A, sclssq );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sclssq->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCORE_zsytf2_nopiv( cham_uplo_t uplo, int n, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}
#endif

int
TCORE_ztplqt( int                    M,
              int                    N,
              int                    L,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          B,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( A, B, T );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztpmlqt( cham_side_t            side,
               cham_trans_t           trans,
               int                    M,
               int                    N,
               int                    K,
               int                    L,
               int                    IB,
               const CHAM_tile_t *    V,
               const CHAM_tile_t *    T,
               CHAM_tile_t *          A,
               CHAM_tile_t *          B,
               CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( V, T, A, B );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztpmqrt( cham_side_t            side,
               cham_trans_t           trans,
               int                    M,
               int                    N,
               int                    K,
               int                    L,
               int                    IB,
               const CHAM_tile_t *    V,
               const CHAM_tile_t *    T,
               CHAM_tile_t *          A,
               CHAM_tile_t *          B,
               CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( V, T, A, B );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztpqrt( int                    M,
              int                    N,
              int                    L,
              int                    IB,
              CHAM_tile_t *          A,
              CHAM_tile_t *          B,
              CHAM_tile_t *          T,
              CHAMELEON_Complex64_t *WORK )
{
    coreblas_kernel_trace( A, B, T );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztradd( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   M,
              int                   N,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         B )
{
    coreblas_kernel_trace( A, B );
    return 0;
}

void
TCORE_ztrasm( __attribute__((unused)) cham_store_t       storev,
              __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) cham_diag_t        diag,
              __attribute__((unused)) int                M,
              __attribute__((unused)) int                N,
              const CHAM_tile_t *A,
              __attribute__((unused)) double *           work )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_ztrmm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAM_tile_t *         B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

void
TCORE_ztrsm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   M,
             int                   N,
             CHAMELEON_Complex64_t alpha,
             const CHAM_tile_t *   A,
             CHAM_tile_t *         B )
{
    coreblas_kernel_trace( A, B );
}

int
TCORE_ztrssq( cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t *      sclssq )
{
    coreblas_kernel_trace( A, sclssq );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sclssq->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_ztrtri( cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t *A, int *info )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
}

int
TCORE_ztsmlq_hetra1( cham_side_t            side,
                     cham_trans_t           trans,
                     int                    m1,
                     int                    n1,
                     int                    m2,
                     int                    n2,
                     int                    k,
                     int                    ib,
                     CHAM_tile_t *          A1,
                     CHAM_tile_t *          A2,
                     const CHAM_tile_t *    V,
                     const CHAM_tile_t *    T,
                     CHAMELEON_Complex64_t *WORK,
                     int                    ldwork )
{
    coreblas_kernel_trace( A1, A2, V, T );
    assert( A1->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A2->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztsmqr_hetra1( cham_side_t            side,
                     cham_trans_t           trans,
                     int                    m1,
                     int                    n1,
                     int                    m2,
                     int                    n2,
                     int                    k,
                     int                    ib,
                     CHAM_tile_t *          A1,
                     CHAM_tile_t *          A2,
                     const CHAM_tile_t *    V,
                     const CHAM_tile_t *    T,
                     CHAMELEON_Complex64_t *WORK,
                     int                    ldwork )
{
    coreblas_kernel_trace( A1, A2, V, T );
    assert( A1->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A2->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_ztstrf( int                    M,
              int                    N,
              int                    IB,
              int                    NB,
              CHAM_tile_t *          U,
              CHAM_tile_t *          A,
              CHAM_tile_t *          L,
              int *                  IPIV,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK,
              int *                  INFO )
{
    coreblas_kernel_trace( U, A, L );
    assert( U->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( L->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zunmlq( cham_side_t            side,
              cham_trans_t           trans,
              int                    M,
              int                    N,
              int                    K,
              int                    IB,
              const CHAM_tile_t *    V,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK )
{
    coreblas_kernel_trace( V, T, C );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zunmqr( cham_side_t            side,
              cham_trans_t           trans,
              int                    M,
              int                    N,
              int                    K,
              int                    IB,
              const CHAM_tile_t *    V,
              const CHAM_tile_t *    T,
              CHAM_tile_t *          C,
              CHAMELEON_Complex64_t *WORK,
              int                    LDWORK )
{
    coreblas_kernel_trace( V, T, C );
    assert( V->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( T->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgesum( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sum )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sum->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zcesca( int center,
              int scale,
              cham_store_t axis,
              int                M,
              int                N,
              int                Mt,
              int                Nt,
              const CHAM_tile_t *Gi,
              const CHAM_tile_t *Gj,
              const CHAM_tile_t *G,
              const CHAM_tile_t *Di,
              const CHAM_tile_t *Dj,
              CHAM_tile_t *      A )
{
    assert( Gi->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( Gj->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( G->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( Di->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( Dj->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

int
TCORE_zgram( cham_uplo_t        uplo,
             int                M,
             int                N,
             int                Mt,
             int                Nt,
             const CHAM_tile_t *Di,
             const CHAM_tile_t *Dj,
             const CHAM_tile_t *D,
             CHAM_tile_t *      A )
{
    coreblas_kernel_trace( Di, Dj, D, A );
    assert( Di->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( Dj->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( D->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return 0;
}

void
TCORE_zprint( __attribute__((unused)) FILE *file,
              __attribute__((unused)) const char *header,
              __attribute__((unused)) cham_uplo_t uplo,
              __attribute__((unused)) int M,
              __attribute__((unused)) int N,
              __attribute__((unused)) int Am,
              __attribute__((unused)) int An,
              const CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
}
