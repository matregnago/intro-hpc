/**
 *
 * @file cuda_htile_empty.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA additional precision kernel empty interface from
 *        CHAM_tile_t layout.
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#include "gpucublas.h"

#if defined(GPUCUBLAS_HAVE_CUBLASHGEMM)
int
TCUDA_hgemm( __attribute__((unused)) cham_trans_t transa,
             __attribute__((unused)) cham_trans_t transb,
             __attribute__((unused)) int m,
             __attribute__((unused)) int n,
             __attribute__((unused)) int k,
             __attribute__((unused)) const CHAMELEON_Real16_t *alpha,
             const CHAM_tile_t *A,
             const CHAM_tile_t *B,
             __attribute__((unused)) const CHAMELEON_Real16_t *beta,
             CHAM_tile_t *C,
             __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}
#endif

#if defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
int
TCUDA_gemmex( __attribute__((unused)) cham_trans_t transa,
              __attribute__((unused)) cham_trans_t transb,
              __attribute__((unused)) int m,
              __attribute__((unused)) int n,
              __attribute__((unused)) int k,
              __attribute__((unused)) const void *alpha,
              const CHAM_tile_t *A,
              const CHAM_tile_t *B,
              __attribute__((unused)) const void *beta,
              CHAM_tile_t *C,
              __attribute__((unused)) cublasHandle_t handle )
{
    gpucublas_kernel_trace( A, B, C );
    return 0;
}
#endif
