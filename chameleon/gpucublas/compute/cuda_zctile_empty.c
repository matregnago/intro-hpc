/**
 *
 * @file cuda_zctile_empty.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA mixed-precision kernel empty interface from CHAM_tile_t layout.
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
TCUDA_clag2z( __attribute__((unused)) int m,
              __attribute__((unused)) int n,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}

int
TCUDA_zlag2c( __attribute__((unused)) int m,
              __attribute__((unused)) int n,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B,
              __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B );
    return 0;
}
