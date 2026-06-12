/**
 *
 * @file cuda_ztile_empty.c
 *
 * @copyright 2025-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA kernel empty interface from CHAM_tile_t layout.
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
TCUDA_zgeadd( __attribute__((unused)) cham_trans_t           trans,
              __attribute__((unused)) int                    m,
              __attribute__((unused)) int                    n,
              __attribute__((unused)) const cuDoubleComplex *alpha,
              const CHAM_tile_t *A,
              __attribute__((unused)) const cuDoubleComplex *beta,
              CHAM_tile_t       *B,
              __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zgemerge( __attribute__((unused)) cham_side_t     side,
                __attribute__((unused)) cham_diag_t     diag,
                __attribute__((unused)) int             M,
                __attribute__((unused)) int             N,
                const CHAM_tile_t                      *A,
                CHAM_tile_t                            *B,
                __attribute__((unused)) cublasHandle_t  handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zgemm( __attribute__((unused)) cham_trans_t           transA,
             __attribute__((unused)) cham_trans_t           transB,
             __attribute__((unused)) int                    m,
             __attribute__((unused)) int                    n,
             __attribute__((unused)) int                    k,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t                             *A,
             const CHAM_tile_t                             *B,
             __attribute__((unused)) const cuDoubleComplex *beta,
             CHAM_tile_t       *C,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}

#if defined( PRECISION_z ) || defined( PRECISION_c )
int
TCUDA_zhemm( __attribute__((unused)) cham_side_t           side,
             __attribute__((unused)) cham_uplo_t           uplo,
             __attribute__((unused)) int                   m,
             __attribute__((unused)) int                   n,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t *A,
             const CHAM_tile_t *B,
             __attribute__((unused)) const cuDoubleComplex *beta,
             CHAM_tile_t       *C,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}

int
TCUDA_zher2k( __attribute__((unused)) cham_uplo_t           uplo,
              __attribute__((unused)) cham_trans_t          trans,
              __attribute__((unused)) int                   n,
              __attribute__((unused)) int                   k,
              __attribute__((unused)) const cuDoubleComplex *alpha,
              const CHAM_tile_t *A,
              const CHAM_tile_t *B,
              __attribute__((unused)) const double          *beta,
              CHAM_tile_t       *C,
              __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}

int
TCUDA_zherk( __attribute__((unused)) cham_uplo_t       uplo,
             __attribute__((unused)) cham_trans_t      trans,
             __attribute__((unused)) int               n,
             __attribute__((unused)) int               k,
             __attribute__((unused)) const double      *alpha,
             const CHAM_tile_t *A,
             __attribute__((unused)) const double      *beta,
             CHAM_tile_t       *C,
             __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, C );
    return 0;
}
#endif

int
TCUDA_zherfb( __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) int                n,
              __attribute__((unused)) int                k,
              __attribute__((unused)) int                ib,
              __attribute__((unused)) int                nb,
              const CHAM_tile_t *A,
              const CHAM_tile_t *T,
              CHAM_tile_t *      C,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                ldwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    return 0;
}

int
TCUDA_zlarfb( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) cham_dir_t         direct,
              __attribute__((unused)) cham_store_t       storev,
              __attribute__((unused)) int                M,
              __attribute__((unused)) int                N,
              __attribute__((unused)) int                K,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      C,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                ldwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, C, WORK );
    return 0;
}

int
TCUDA_zlacpy( __attribute__((unused)) cham_uplo_t    uplo,
              __attribute__((unused)) int            M,
              __attribute__((unused)) int            N,
              const CHAM_tile_t                     *A,
              CHAM_tile_t                           *B,
              __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zlacpyx( __attribute__((unused)) cham_uplo_t    uplo,
               __attribute__((unused)) int            M,
               __attribute__((unused)) int            N,
               __attribute__((unused)) int            displA,
               const CHAM_tile_t                     *A,
               __attribute__((unused)) int            LDA,
               __attribute__((unused)) int            displB,
               CHAM_tile_t                           *B,
               __attribute__((unused)) int            LDB,
               __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zlaset( __attribute__((unused)) cham_uplo_t            uplo,
              __attribute__((unused)) int                    m,
              __attribute__((unused)) int                    n,
              __attribute__((unused)) const cuDoubleComplex *alpha,
              __attribute__((unused)) const cuDoubleComplex *beta,
              CHAM_tile_t                                   *A,
              __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A );
    return 0;
}

int
TCUDA_zlatro( __attribute__((unused)) cham_uplo_t    uplo,
              __attribute__((unused)) cham_trans_t   trans,
              __attribute__((unused)) int            M,
              __attribute__((unused)) int            N,
              CHAM_tile_t                           *A,
              CHAM_tile_t                           *B,
              __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zparfb( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) cham_dir_t         direct,
              __attribute__((unused)) cham_store_t       storev,
              __attribute__((unused)) int                M1,
              __attribute__((unused)) int                N1,
              __attribute__((unused)) int                M2,
              __attribute__((unused)) int                N2,
              __attribute__((unused)) int                K,
              __attribute__((unused)) int                L,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    return 0;
}

int
TCUDA_zpotrf( __attribute__((unused)) cham_uplo_t        uplo,
              __attribute__((unused)) int                n,
              CHAM_tile_t                               *A,
              cuComplexDouble                           *WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) int               *d_info,
              __attribute__((unused)) cusolverDnHandle_t handle )
{
    gpucublas_kernel_trace( A, WORK );
    return 0;
}

int
TCUDA_zsymm( __attribute__((unused)) cham_side_t           side,
             __attribute__((unused)) cham_uplo_t           uplo,
             __attribute__((unused)) int                   m,
             __attribute__((unused)) int                   n,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t *A,
             const CHAM_tile_t *B,
             __attribute__((unused)) const cuDoubleComplex *beta,
             CHAM_tile_t       *C,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}

int
TCUDA_zsyr2k( __attribute__((unused)) cham_uplo_t           uplo,
              __attribute__((unused)) cham_trans_t          trans,
              __attribute__((unused)) int                   n,
              __attribute__((unused)) int                   k,
              __attribute__((unused)) const cuDoubleComplex *alpha,
              const CHAM_tile_t *A,
              const CHAM_tile_t *B,
              __attribute__((unused)) const cuDoubleComplex *beta,
              CHAM_tile_t       *C,
              __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}

int
TCUDA_zsyrk( __attribute__((unused)) cham_uplo_t           uplo,
             __attribute__((unused)) cham_trans_t          trans,
             __attribute__((unused)) int                   n,
             __attribute__((unused)) int                   k,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t *A,
             __attribute__((unused)) const cuDoubleComplex *beta,
             CHAM_tile_t       *C,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, C );
    return 0;
}

int
TCUDA_ztpmqrt( __attribute__((unused)) cham_side_t        side,
               __attribute__((unused)) cham_trans_t       trans,
               __attribute__((unused)) int                M,
               __attribute__((unused)) int                N,
               __attribute__((unused)) int                K,
               __attribute__((unused)) int                L,
               __attribute__((unused)) int                IB,
               const CHAM_tile_t *V,
               const CHAM_tile_t *T,
               CHAM_tile_t       *A,
               CHAM_tile_t       *B,
               CHAM_tile_t       *WORK,
               __attribute__((unused)) int                lwork,
               __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, A, B, WORK );
    return 0;
}

int
TCUDA_ztpmlqt( __attribute__((unused)) cham_side_t        side,
               __attribute__((unused)) cham_trans_t       trans,
               __attribute__((unused)) int                M,
               __attribute__((unused)) int                N,
               __attribute__((unused)) int                K,
               __attribute__((unused)) int                L,
               __attribute__((unused)) int                IB,
               const CHAM_tile_t *V,
               const CHAM_tile_t *T,
               CHAM_tile_t       *A,
               CHAM_tile_t       *B,
               CHAM_tile_t       *WORK,
               __attribute__((unused)) int                lwork,
               __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( V, T, A, B, WORK );
    return 0;
}

int
TCUDA_ztrmm( __attribute__((unused)) cham_side_t           side,
             __attribute__((unused)) cham_uplo_t           uplo,
             __attribute__((unused)) cham_trans_t          transA,
             __attribute__((unused)) cham_diag_t           diag,
             __attribute__((unused)) int                   m,
             __attribute__((unused)) int                   n,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t *A,
             CHAM_tile_t       *B,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_ztrsm( __attribute__((unused)) cham_side_t           side,
             __attribute__((unused)) cham_uplo_t           uplo,
             __attribute__((unused)) cham_trans_t          transA,
             __attribute__((unused)) cham_diag_t           diag,
             __attribute__((unused)) int                   m,
             __attribute__((unused)) int                   n,
             __attribute__((unused)) const cuDoubleComplex *alpha,
             const CHAM_tile_t *A,
             CHAM_tile_t       *B,
             __attribute__((unused)) cublasHandle_t         handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_ztsmlq( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) int                M1,
              __attribute__((unused)) int                N1,
              __attribute__((unused)) int                M2,
              __attribute__((unused)) int                N2,
              __attribute__((unused)) int                K,
              __attribute__((unused)) int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    return 0;
}

int
TCUDA_zttmlq( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) int                M1,
              __attribute__((unused)) int                N1,
              __attribute__((unused)) int                M2,
              __attribute__((unused)) int                N2,
              __attribute__((unused)) int                K,
              __attribute__((unused)) int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    return 0;
}

int
TCUDA_ztsmqr( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) int                M1,
              __attribute__((unused)) int                N1,
              __attribute__((unused)) int                M2,
              __attribute__((unused)) int                N2,
              __attribute__((unused)) int                K,
              __attribute__((unused)) int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    return 0;
}

int
TCUDA_zttmqr( __attribute__((unused)) cham_side_t        side,
              __attribute__((unused)) cham_trans_t       trans,
              __attribute__((unused)) int                M1,
              __attribute__((unused)) int                N1,
              __attribute__((unused)) int                M2,
              __attribute__((unused)) int                N2,
              __attribute__((unused)) int                K,
              __attribute__((unused)) int                IB,
              CHAM_tile_t *      A1,
              CHAM_tile_t *      A2,
              const CHAM_tile_t *V,
              const CHAM_tile_t *T,
              CHAM_tile_t *      WORK,
              __attribute__((unused)) int                lwork,
              __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A1, A2, V, T, WORK );
    return 0;
}

int
TCUDA_zunmlqt( __attribute__((unused)) cham_side_t        side,
               __attribute__((unused)) cham_trans_t       trans,
               __attribute__((unused)) int                M,
               __attribute__((unused)) int                N,
               __attribute__((unused)) int                K,
               __attribute__((unused)) int                IB,
               const CHAM_tile_t *A,
               const CHAM_tile_t *T,
               CHAM_tile_t *      C,
               CHAM_tile_t *      WORK,
               __attribute__((unused)) int                ldwork,
               __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    return 0;
}

int
TCUDA_zunmqrt( __attribute__((unused)) cham_side_t        side,
               __attribute__((unused)) cham_trans_t       trans,
               __attribute__((unused)) int                M,
               __attribute__((unused)) int                N,
               __attribute__((unused)) int                K,
               __attribute__((unused)) int                IB,
               const CHAM_tile_t *A,
               const CHAM_tile_t *T,
               CHAM_tile_t *      C,
               CHAM_tile_t *      WORK,
               __attribute__((unused)) int                ldwork,
               __attribute__((unused)) cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, T, C, WORK );
    return 0;
}
