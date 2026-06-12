/**
 *
 * @file core_zctile_empty.c
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
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#include "coreblas.h"
#include "coreblas/coreblas_zctile.h"

void
TCORE_clag2z( int M, int N,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
}

void
TCORE_zlag2c( int M, int N,
              const CHAM_tile_t *A,
              CHAM_tile_t       *B, int *info )
{
    coreblas_kernel_trace( A, B );
    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );
}
