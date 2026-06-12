/**
 *
 * @file cuda_ztile.c
 *
 * @copyright 2025-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"
#include "gpucublas/gpucublas_ztile.h"

int
TCUDA_zgeadd( cham_trans_t           trans,
              int                    m,
              int                    n,
              const cuDoubleComplex *alpha,
              const CHAM_tile_t     *A,
              const cuDoubleComplex *beta,
              CHAM_tile_t           *B,
              cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zgeadd( trans, m, n, alpha, (const cuDoubleComplex *)A->mat, A->ld, beta, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_zgemerge( cham_side_t    side,
                cham_diag_t    diag,
                int            M,
                int            N,
                const CHAM_tile_t *A,
                CHAM_tile_t *      B,
                cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zgemerge( side, diag, M, N, (const cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_zgemm( cham_trans_t           transA,
             cham_trans_t           transB,
             int                    m,
             int                    n,
             int                    k,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             const CHAM_tile_t     *B,
             const cuDoubleComplex *beta,
             CHAM_tile_t           *C,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zgemm( transA, transB, m, n, k, alpha, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)B->mat, B->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCUDA_zhemm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   m,
             int                   n,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             const CHAM_tile_t     *B,
             const cuDoubleComplex *beta,
             CHAM_tile_t           *C,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zhemm( side, uplo, m, n, alpha, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)B->mat, B->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

int
TCUDA_zher2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   n,
              int                   k,
              const cuDoubleComplex *alpha,
              const CHAM_tile_t     *A,
              const CHAM_tile_t     *B,
              const double          *beta,
              CHAM_tile_t           *C,
              cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zher2k( uplo, trans, n, k, alpha, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)B->mat, B->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

int
TCUDA_zherk( cham_uplo_t       uplo,
             cham_trans_t      trans,
             int               n,
             int               k,
             const double      *alpha,
             const CHAM_tile_t *A,
             const double      *beta,
             CHAM_tile_t       *C,
             cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zherk( uplo, trans, n, k, alpha, (const cuDoubleComplex *)A->mat, A->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}
#endif

int
TCUDA_zherfb( cham_uplo_t        uplo,
              int                n,
              int                k,
              int                ib,
              int                nb,
              const CHAM_tile_t *A,
              const CHAM_tile_t *T,
              CHAM_tile_t *      C,
              CHAM_tile_t *      WORK,
              int                ldwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zherfb( uplo, n, k, ib, nb, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)C->mat, C->ld, (cuDoubleComplex *)WORK->mat, ldwork, handle );
}

int
TCUDA_zlarfb( cham_side_t        side,
              cham_trans_t       trans,
              cham_dir_t         direct,
              cham_store_t       storev,
              int                M,
              int                N,
              int                K,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      C,
              CHAM_tile_t *      WORK,
              int                ldwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, C, WORK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zlarfb( side, trans, direct, storev, M, N, K, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)C->mat, C->ld, (cuDoubleComplex *)WORK->mat, ldwork, handle );
}

int
TCUDA_zlacpy( cham_uplo_t        uplo,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zlacpy( uplo, M, N, (const cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_zlacpyx( cham_uplo_t        uplo,
               int                M,
               int                N,
               int                displA,
               const CHAM_tile_t *A,
               int                LDA,
               int                displB,
               CHAM_tile_t       *B,
               int                LDB,
               cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );

    const cuDoubleComplex *Aptr = (const cuDoubleComplex *)A->mat;
    cuDoubleComplex *Bptr = (cuDoubleComplex *)B->mat;
    return CUDA_zlacpy( uplo, M, N, Aptr + displA, LDA, Bptr + displB, LDB, handle );
}

int
TCUDA_zlaset( cham_uplo_t            uplo,
              int                    m,
              int                    n,
              const cuDoubleComplex *alpha,
              const cuDoubleComplex *beta,
              CHAM_tile_t           *A,
              cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zlaset( uplo, m, n, alpha, beta, CHAM_tile_get_ptr( A ), A->ld, handle );
}

int
TCUDA_zlatro( cham_uplo_t        uplo,
              cham_trans_t       trans,
              int                M,
              int                N,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zlatro( uplo, trans, M, N, (cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_zparfb( cham_side_t        side,
              cham_trans_t       trans,
              cham_dir_t         direct,
              cham_store_t       storev,
              int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                L,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              int                lwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zparfb( side, trans, direct, storev, M1, N1, M2, N2, K, L, (cuDoubleComplex *)A1->mat, A1->ld, (cuDoubleComplex *)A2->mat, A2->ld, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_zpotrf( cham_uplo_t        uplo,
              int                n,
              CHAM_tile_t       *A,
              cuDoubleComplex   *dW,
              int                lwork,
              int               *d_info,
              cusolverDnHandle_t handle )
{
    gpucublas_kernel_trace( A );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zpotrf( uplo, n, (cuDoubleComplex *)A->mat, A->ld, dW, lwork, d_info, handle );
}

int
TCUDA_zsymm( cham_side_t           side,
             cham_uplo_t           uplo,
             int                   m,
             int                   n,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             const CHAM_tile_t     *B,
             const cuDoubleComplex *beta,
             CHAM_tile_t           *C,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zsymm( side, uplo, m, n, alpha, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)B->mat, B->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

int
TCUDA_zsyr2k( cham_uplo_t           uplo,
              cham_trans_t          trans,
              int                   n,
              int                   k,
              const cuDoubleComplex *alpha,
              const CHAM_tile_t     *A,
              const CHAM_tile_t     *B,
              const cuDoubleComplex *beta,
              CHAM_tile_t           *C,
              cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zsyr2k( uplo, trans, n, k, alpha, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)B->mat, B->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

int
TCUDA_zsyrk( cham_uplo_t           uplo,
             cham_trans_t          trans,
             int                   n,
             int                   k,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             const cuDoubleComplex *beta,
             CHAM_tile_t           *C,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zsyrk( uplo, trans, n, k, alpha, (const cuDoubleComplex *)A->mat, A->ld, beta, (cuDoubleComplex *)C->mat, C->ld, handle );
}

int
TCUDA_ztpmqrt( cham_side_t        side,
               cham_trans_t       trans,
               int                M,
               int                N,
               int                K,
               int                L,
               int                IB,
               const CHAM_tile_t *V,
               const CHAM_tile_t *T,
               CHAM_tile_t       *A,
               CHAM_tile_t       *B,
               CHAM_tile_t       *WORK,
               int                lwork,
               cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, A, B, WORK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztpmqrt( side, trans, M, N, K, L, IB, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_ztpmlqt( cham_side_t        side,
               cham_trans_t       trans,
               int                M,
               int                N,
               int                K,
               int                L,
               int                IB,
               const CHAM_tile_t *V,
               const CHAM_tile_t *T,
               CHAM_tile_t       *A,
               CHAM_tile_t       *B,
               CHAM_tile_t       *WORK,
               int                lwork,
               cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, A, B, WORK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztpmlqt( side, trans, M, N, K, L, IB, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_ztrmm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   m,
             int                   n,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             CHAM_tile_t           *B,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztrmm( side, uplo, transA, diag, m, n, alpha, (const cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_ztrsm( cham_side_t           side,
             cham_uplo_t           uplo,
             cham_trans_t          transA,
             cham_diag_t           diag,
             int                   m,
             int                   n,
             const cuDoubleComplex *alpha,
             const CHAM_tile_t     *A,
             CHAM_tile_t           *B,
             cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztrsm( side, uplo, transA, diag, m, n, alpha, (const cuDoubleComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_ztsmlq( cham_side_t        side,
              cham_trans_t       trans,
              int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              int                lwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztsmlq( side, trans, M1, N1, M2, N2, K, IB, (cuDoubleComplex *)A1->mat, A1->ld, (cuDoubleComplex *)A2->mat, A2->ld, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_zttmlq( cham_side_t        side,
              cham_trans_t       trans,
              int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              int                lwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zttmlq( side, trans, M1, N1, M2, N2, K, IB, (cuDoubleComplex *)A1->mat, A1->ld, (cuDoubleComplex *)A2->mat, A2->ld, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_ztsmqr( cham_side_t        side,
              cham_trans_t       trans,
              int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              int                lwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_ztsmqr( side, trans, M1, N1, M2, N2, K, IB, (cuDoubleComplex *)A1->mat, A1->ld, (cuDoubleComplex *)A2->mat, A2->ld, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_zttmqr( cham_side_t        side,
              cham_trans_t       trans,
              int                M1,
              int                N1,
              int                M2,
              int                N2,
              int                K,
              int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              int                lwork,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    assert( A1->format & CHAMELEON_TILE_FULLRANK );
    assert( A2->format & CHAMELEON_TILE_FULLRANK );
    assert( V->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zttmqr( side, trans, M1, N1, M2, N2, K, IB, (cuDoubleComplex *)A1->mat, A1->ld, (cuDoubleComplex *)A2->mat, A2->ld, (const cuDoubleComplex *)V->mat, V->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)WORK->mat, lwork, handle );
}

int
TCUDA_zunmlqt( cham_side_t        side,
               cham_trans_t       trans,
               int                M,
               int                N,
               int                K,
               int                IB,
               const CHAM_tile_t *A,
               const CHAM_tile_t *T,
               CHAM_tile_t *      C,
               CHAM_tile_t *      WORK,
               int                ldwork,
               cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zunmlqt( side, trans, M, N, K, IB, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)C->mat, C->ld, (cuDoubleComplex *)WORK->mat, ldwork, handle );
}

int
TCUDA_zunmqrt( cham_side_t        side,
               cham_trans_t       trans,
               int                M,
               int                N,
               int                K,
               int                IB,
               const CHAM_tile_t *A,
               const CHAM_tile_t *T,
               CHAM_tile_t *      C,
               CHAM_tile_t *      WORK,
               int                ldwork,
               cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( T->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    assert( WORK->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zunmqrt( side, trans, M, N, K, IB, (const cuDoubleComplex *)A->mat, A->ld, (const cuDoubleComplex *)T->mat, T->ld, (cuDoubleComplex *)C->mat, C->ld, (cuDoubleComplex *)WORK->mat, ldwork, handle );
}
