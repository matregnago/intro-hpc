/**
 *
 * @file cuda_htile.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA additional precision kernel interface from CHAM_tile_t
 *        layout to the real one.
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#include "gpucublas.h"

#if defined(GPUCUBLAS_HAVE_CUBLASHGEMM)
int
TCUDA_hgemm( cham_trans_t transa, cham_trans_t transb,
             int m, int n, int k,
             const CHAMELEON_Real16_t *alpha,
             const CHAM_tile_t *A,
             const CHAM_tile_t *B,
             const CHAMELEON_Real16_t *beta,
             CHAM_tile_t *C,
             cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_hgemm( transa, transb, m, n, k, alpha, A->mat, A->ld, B->mat, B->ld, beta, C->mat, C->ld, handle );
}
#endif

#if defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
int
TCUDA_gemmex( cham_trans_t transa, cham_trans_t transb,
              int m, int n, int k,
              const void *alpha,
              const CHAM_tile_t *A,
              const CHAM_tile_t *B,
              const void *beta,
              CHAM_tile_t *C,
              cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B, C );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    assert( C->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_gemmex( transa, transb, m, n, k, alpha, A->mat, A->ld, A->flttype, B->mat, B->ld, B->flttype, beta, C->mat, C->ld, C->flttype, handle );
}
#endif
