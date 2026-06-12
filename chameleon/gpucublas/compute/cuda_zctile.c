/**
 *
 * @file cuda_zctile.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA mixed-precision kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions mixed zc -> ds
 *
 */
#include "gpucublas.h"
#include "gpucublas/gpucublas_zctile.h"

int
TCUDA_clag2z( int m, int n,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_clag2z( m, n, (const cuFloatComplex *)A->mat, A->ld, (cuDoubleComplex *)B->mat, B->ld, handle );
}

int
TCUDA_zlag2c( int m, int n,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              cublasHandle_t     handle )
{
    gpucublas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
    return CUDA_zlag2c( m, n, (const cuDoubleComplex *)A->mat, A->ld, (cuFloatComplex *)B->mat, B->ld, handle );
}
