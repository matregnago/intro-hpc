/**
 *
 * @file core_ztile.c
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
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Abel Calluaud
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include "coreblas.h"
#include "coreblas/coreblas_ztile.h"

#if defined( CHAMELEON_USE_HMATOSS )
#include "coreblas/hmat.h"
#endif

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_dlag2z( cham_uplo_t uplo, int M, int N,
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
    return CORE_zaxpy( M, alpha, CHAM_tile_get_ptr( A ), incA, CHAM_tile_get_ptr( B ), incB );
}

int
TCORE_zgeadd( cham_trans_t          trans,
              int                   M,
              int                   N,
              CHAMELEON_Complex64_t alpha,
              const CHAM_tile_t *   A,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         B )
{
    coreblas_kernel_trace( A, B );
    if ( (A->format & CHAMELEON_TILE_DESC) &&
         (B->format & CHAMELEON_TILE_DESC) )
    {
        assert(0);
    }

    return CORE_zgeadd( trans, M, N,
                        alpha, CHAM_tile_get_ptr( A ), A->ld,
                        beta,  CHAM_tile_get_ptr( B ), B->ld );
}

int
TCORE_zgelqt( int                    M,
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
    return CORE_zgelqt( M, N, IB, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( T ), T->ld, TAU, WORK );
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
    CORE_zgemv(
        trans, M, N, alpha, A->mat, A->ld, x->mat, incX, beta, y->mat, incY );
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
    if ( ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) &&
         ( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) &&
         ( C->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) )
    {
        CORE_zgemm( transA, transB, M, N, K, alpha,
                    CHAM_tile_get_ptr( A ), A->ld,
                    CHAM_tile_get_ptr( B ), B->ld,
                    beta,
                    CHAM_tile_get_ptr( C ), C->ld );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( ( A->format & B->format & C->format & CHAMELEON_TILE_HMAT ) )
    {
        hmat_zgemm( chameleon_lapack_const( transA ),
                    chameleon_lapack_const( transB ),
                    &alpha, A->mat, B->mat,
                    &beta, C->mat );
    }
    else if ( ( A->format & CHAMELEON_TILE_HMAT     ) &&
              ( B->format & C->format & CHAMELEON_TILE_FULLRANK ) )
    {
        assert( transB == ChamNoTrans );
        hmat_zgemv( chameleon_lapack_const( transA ),
                    &alpha, A->mat,
                    CHAM_tile_get_ptr( B ),
                    &beta,
                    CHAM_tile_get_ptr( C ), C->n );
    }
#endif
    else {
        assert( 0 );
    }
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
    return CORE_zgeqrt( M, N, IB, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( T ), T->ld, TAU, WORK );
}

int
TCORE_zgessm( int M, int N, int K, int IB, const int *IPIV, const CHAM_tile_t *L, CHAM_tile_t *A )
{
    coreblas_kernel_trace( L, A );
    assert( L->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zgessm( M, N, K, IB, IPIV, CHAM_tile_get_ptr( L ), L->ld, CHAM_tile_get_ptr( A ), A->ld );
}

int
TCORE_zgessq( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sclssq )
{
    coreblas_kernel_trace( A, sclssq );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sclssq->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zgessq( storev, M, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( sclssq ) );
}

int
TCORE_zgetrf( int M, int N, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    coreblas_kernel_trace( A );

    int rc = -1;
    if ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) {
        rc = CORE_zgetrf( M, N, CHAM_tile_get_ptr( A ), A->ld, IPIV, INFO );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( A->format & CHAMELEON_TILE_HMAT ) {
        rc = hmat_zgetrf( A->mat );
        assert( rc == 0 );
    }
#endif
    else {
        assert( 0 );
    }
    return rc;
}

int
TCORE_zgetrf_incpiv( int M, int N, int IB, CHAM_tile_t *A, int *IPIV, int *INFO )
{
    coreblas_kernel_trace( A );

    if ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) {
        return CORE_zgetrf_incpiv( M, N, IB, CHAM_tile_get_ptr( A ), A->ld, IPIV, INFO );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( A->format & CHAMELEON_TILE_HMAT ) {
        return hmat_zgetrf( A->mat );
    }
#endif
    else {
        assert( 0 );
    }
    return -1;
}

int
TCORE_zgetrf_nopiv( int M, int N, int IB, CHAM_tile_t *A, int *INFO )
{
    coreblas_kernel_trace( A );

    int rc = -1;
    *INFO  = 0;

    if ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) {
        rc = CORE_zgetrf_nopiv( M, N, IB, CHAM_tile_get_ptr( A ), A->ld, INFO );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( A->format & CHAMELEON_TILE_HMAT ) {
        rc     = hmat_zgetrf( A->mat );
        assert( rc == 0 );
    }
#endif
    else {
        assert( 0 );
    }
    return rc;
}

void
TCORE_zhe2ge( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zhe2ge( uplo, M, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld );
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
    CORE_zhemm( side, uplo, M, N, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    CORE_zherk( uplo, trans, N, K, alpha, CHAM_tile_get_ptr( A ), A->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    CORE_zher2k( uplo, trans, N, K, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    return CORE_zherfb(
        uplo, N, K, IB, NB, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( T ), T->ld, CHAM_tile_get_ptr( C ), C->ld, WORK, ldwork );
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
    return CORE_zhessq( storev, uplo, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( sclssq ) );
}
#endif

void
TCORE_zlacpy( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B )
{
    if (( A->format & CHAMELEON_TILE_DESC ) &&
        ( B->format & CHAMELEON_TILE_DESC ) )
    {
        assert(0); /* This should have been handled at the codelet level */
    }
    CORE_zlacpy( uplo, M, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld );
}

void
TCORE_zlacpyx( cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t *A, int LDA, int displB, CHAM_tile_t *B, int LDB )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );

    const CHAMELEON_Complex64_t *Aptr = CHAM_tile_get_ptr( A );
    CHAMELEON_Complex64_t       *Bptr = CHAM_tile_get_ptr( B );
    CORE_zlacpy( uplo, M, N, Aptr + displA, LDA, Bptr + displB, LDB );
}

void
TCORE_zlange( cham_normtype_t    norm,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlange( norm, M, N, CHAM_tile_get_ptr( A ), A->ld, work, normA );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
void
TCORE_zlanhe( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlanhe( norm, uplo, N, CHAM_tile_get_ptr( A ), A->ld, work, normA );
}
#endif

void
TCORE_zlansy( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlansy( norm, uplo, N, CHAM_tile_get_ptr( A ), A->ld, work, normA );
}

void
TCORE_zlantr( cham_normtype_t    norm,
              cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work,
              double *           normA )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlantr( norm, uplo, diag, M, N, CHAM_tile_get_ptr( A ), A->ld, work, normA );
}

int
TCORE_zlascal( cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zlascal( uplo, m, n, alpha, CHAM_tile_get_ptr( A ), A->ld );
}

void
TCORE_zlaset( cham_uplo_t           uplo,
              int                   m,
              int                   n,
              CHAMELEON_Complex64_t alpha,
              CHAMELEON_Complex64_t beta,
              CHAM_tile_t *         A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlaset( uplo, m, n, alpha, beta, CHAM_tile_get_ptr( A ), A->ld );
}

void
TCORE_zlaset2( cham_uplo_t uplo, int m, int n, CHAMELEON_Complex64_t alpha, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlaset2( uplo, m, n, alpha, CHAM_tile_get_ptr( A ), A->ld );
}

int
TCORE_zlaswp_get( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *perm )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zlaswp_get( side, m0, m, n, k, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, perm );
}

int
TCORE_zlaswp_get_idx( cham_side_t side, int m0, int m, int n, int k,
                      const CHAM_tile_t *A,
                      CHAM_tile_t       *B,
                      int perm_m, int perm_mt, const int *perm_idx )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zlaswp_get_idx( side, m0, m, n, k,
                                CHAM_tile_get_ptr( A ), A->ld,
                                CHAM_tile_get_ptr( B ), B->ld,
                                perm_m, perm_mt, perm_idx );
}

int
TCORE_zlaswp_set( cham_side_t side, int m0, int m, int n, int k, CHAM_tile_t *A, CHAM_tile_t *B, const int *invp )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zlaswp_set( side, m0, m, n, k, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, invp );
}

int
TCORE_zlaswp_set_idx( cham_side_t side, int m0, int m, int n, int k,
                      const CHAM_tile_t *A,
                      CHAM_tile_t       *B,
                      int invp_m, int invp_mt, const int *invp_idx )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zlaswp_set_idx( side, m0, m, n, k,
                                CHAM_tile_get_ptr( A ), A->ld,
                                CHAM_tile_get_ptr( B ), B->ld,
                                invp_m, invp_mt, invp_idx );
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
    return CORE_zlatro( uplo, trans, M, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld );
}

void
TCORE_zlauum( cham_uplo_t uplo, int N, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_zlauum( uplo, N, CHAM_tile_get_ptr( A ), A->ld );
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
    CORE_zplghe( bump, m, n, CHAM_tile_get_ptr( A ), A->ld, bigM, m0, n0, seed );
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
    CORE_zplgsy( bump, m, n, CHAM_tile_get_ptr( A ), A->ld, bigM, m0, n0, seed );
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
    CORE_zplrnt( m, n, CHAM_tile_get_ptr( A ), A->ld, bigM, m0, n0, seed );
}

void
TCORE_zpotrf( cham_uplo_t uplo, int n, CHAM_tile_t *A, int *INFO )
{
    coreblas_kernel_trace( A );
    if ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) {
        CORE_zpotrf( uplo, n, CHAM_tile_get_ptr( A ), A->ld, INFO );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( A->format & CHAMELEON_TILE_HMAT ) {
        assert( uplo == ChamLower );
        *INFO = hmat_zpotrf( A->mat );
    }
#endif
    else {
        assert( 0 );
    }
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
    return CORE_zssssm( M1,
                        N1,
                        M2,
                        N2,
                        K,
                        IB,
                        CHAM_tile_get_ptr( A1 ),
                        A1->ld,
                        CHAM_tile_get_ptr( A2 ),
                        A2->ld,
                        CHAM_tile_get_ptr( L1 ),
                        L1->ld,
                        CHAM_tile_get_ptr( L2 ),
                        L2->ld,
                        IPIV );
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
    CORE_zsymm( side, uplo, M, N, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    CORE_zsyrk( uplo, trans, N, K, alpha, CHAM_tile_get_ptr( A ), A->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    CORE_zsyr2k( uplo, trans, N, K, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, beta, CHAM_tile_get_ptr( C ), C->ld );
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
    return CORE_zsyssq( storev, uplo, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( sclssq ) );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCORE_zsytf2_nopiv( cham_uplo_t uplo, int n, CHAM_tile_t *A )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zsytf2_nopiv( uplo, n, CHAM_tile_get_ptr( A ), A->ld );
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
    return CORE_ztplqt( M, N, L, IB, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, CHAM_tile_get_ptr( T ), T->ld, WORK );
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
    return CORE_ztpmlqt( side,
                         trans,
                         M,
                         N,
                         K,
                         L,
                         IB,
                         CHAM_tile_get_ptr( V ),
                         V->ld,
                         CHAM_tile_get_ptr( T ),
                         T->ld,
                         CHAM_tile_get_ptr( A ),
                         A->ld,
                         CHAM_tile_get_ptr( B ),
                         B->ld,
                         WORK );
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
    return CORE_ztpmqrt( side,
                         trans,
                         M,
                         N,
                         K,
                         L,
                         IB,
                         CHAM_tile_get_ptr( V ),
                         V->ld,
                         CHAM_tile_get_ptr( T ),
                         T->ld,
                         CHAM_tile_get_ptr( A ),
                         A->ld,
                         CHAM_tile_get_ptr( B ),
                         B->ld,
                         WORK );
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
    return CORE_ztpqrt( M, N, L, IB, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld, CHAM_tile_get_ptr( T ), T->ld, WORK );
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
    if (( A->format & CHAMELEON_TILE_DESC ) &&
        ( B->format & CHAMELEON_TILE_DESC ) )
    {
        assert(0); /* This should have been handled at the codelet level */
    }
    return CORE_ztradd( uplo, trans, M, N, alpha, CHAM_tile_get_ptr( A ), A->ld, beta, CHAM_tile_get_ptr( B ), B->ld );
}

void
TCORE_ztrasm( cham_store_t       storev,
              cham_uplo_t        uplo,
              cham_diag_t        diag,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              double *           work )
{
    coreblas_kernel_trace( A );
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_ztrasm( storev, uplo, diag, M, N, CHAM_tile_get_ptr( A ), A->ld, work );
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
    CORE_ztrmm( side, uplo, transA, diag, M, N, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld );
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

    if ( ( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) &&
         ( B->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) ) )
    {
        CORE_ztrsm( side, uplo, transA, diag, M, N, alpha, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( B ), B->ld );
    }
#if defined( CHAMELEON_USE_HMATOSS )
    else if ( A->format & CHAMELEON_TILE_HMAT ) {
        assert( !(B->format & CHAMELEON_TILE_DESC) );
        hmat_ztrsm( chameleon_lapack_const( side ),
                    chameleon_lapack_const( uplo ),
                    chameleon_lapack_const( transA ),
                    chameleon_lapack_const( diag ),
                    M, N, &alpha, A->mat,
                    ( B->format & CHAMELEON_TILE_HMAT ),
                    CHAM_tile_get_ptr( B ) );
    }
#endif
    else {
        assert( 0 );
    }
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
    double *W = CHAM_tile_get_ptr( sclssq );
    return CORE_ztrssq( uplo, diag, M, N, CHAM_tile_get_ptr( A ), A->ld, W, W + 1 );
}

void
TCORE_ztrtri( cham_uplo_t uplo, cham_diag_t diag, int N, CHAM_tile_t *A, int *info )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    CORE_ztrtri( uplo, diag, N, CHAM_tile_get_ptr( A ), A->ld, info );
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
    return CORE_ztsmlq_hetra1( side,
                               trans,
                               m1,
                               n1,
                               m2,
                               n2,
                               k,
                               ib,
                               CHAM_tile_get_ptr( A1 ),
                               A1->ld,
                               CHAM_tile_get_ptr( A2 ),
                               A2->ld,
                               CHAM_tile_get_ptr( V ),
                               V->ld,
                               CHAM_tile_get_ptr( T ),
                               T->ld,
                               WORK,
                               ldwork );
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
    return CORE_ztsmqr_hetra1( side,
                               trans,
                               m1,
                               n1,
                               m2,
                               n2,
                               k,
                               ib,
                               CHAM_tile_get_ptr( A1 ),
                               A1->ld,
                               CHAM_tile_get_ptr( A2 ),
                               A2->ld,
                               CHAM_tile_get_ptr( V ),
                               V->ld,
                               CHAM_tile_get_ptr( T ),
                               T->ld,
                               WORK,
                               ldwork );
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
    return CORE_ztstrf(
        M, N, IB, NB, CHAM_tile_get_ptr( U ), U->ld, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( L ), L->ld, IPIV, WORK, LDWORK, INFO );
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
    return CORE_zunmlq(
        side, trans, M, N, K, IB, CHAM_tile_get_ptr( V ), V->ld, CHAM_tile_get_ptr( T ), T->ld, CHAM_tile_get_ptr( C ), C->ld, WORK, LDWORK );
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
    return CORE_zunmqr(
        side, trans, M, N, K, IB, CHAM_tile_get_ptr( V ), V->ld, CHAM_tile_get_ptr( T ), T->ld, CHAM_tile_get_ptr( C ), C->ld, WORK, LDWORK );
}

int
TCORE_zgesum( cham_store_t storev, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *sum )
{
    assert( A->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    assert( sum->format & (CHAMELEON_TILE_FULLRANK | CHAMELEON_TILE_DESC) );
    return CORE_zgesum( storev, M, N, CHAM_tile_get_ptr( A ), A->ld, CHAM_tile_get_ptr( sum ) );
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
    return CORE_zcesca(
        center, scale, axis, M, N, Mt, Nt, CHAM_tile_get_ptr( Gi ), Gi->ld, CHAM_tile_get_ptr( Gj ), Gj->ld, CHAM_tile_get_ptr( G ), CHAM_tile_get_ptr( Di ), Di->ld, CHAM_tile_get_ptr( Dj ), Dj->ld, CHAM_tile_get_ptr( A ), A->ld );
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
    return CORE_zgram(
        uplo, M, N, Mt, Nt, CHAM_tile_get_ptr( Di ), Di->ld, CHAM_tile_get_ptr( Dj ), Dj->ld, CHAM_tile_get_ptr( D ), CHAM_tile_get_ptr( A ), A->ld );
}

void
TCORE_zprint( FILE *file, const char *header,
              cham_uplo_t uplo, int M, int N,
              int Am, int An, const CHAM_tile_t *A )
{
    char *s;
    int   rc;
    coreblas_kernel_trace( A );
    assert( A->format & CHAMELEON_TILE_FULLRANK );

#if !defined(CHAMELEON_KERNELS_TRACE)
    rc = asprintf( &s, "%s (%2d, %2d)", header, Am, An );
#else
    rc = asprintf( &s, "%s %s", header, A->name );
#endif
    assert( rc != -1 );
    CORE_zprint( file, s, uplo, M, N, CHAM_tile_get_ptr( A ), A->ld );
    free( s );
    (void)rc;
}
