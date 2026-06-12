/**
 *
 * @file gpucublas_zctile.h
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA mixed-precision tile-based kernel interface header
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions mixed zc -> ds
 *
 */
#ifndef _gpucublas_zctile_h_
#define _gpucublas_zctile_h_

#include "chameleon/config.h"
#include "chameleon/struct.h"

int TCUDA_clag2z( int m, int n, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_zlag2c( int m, int n, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );

#endif /* _gpucublas_zctile_h_ */
