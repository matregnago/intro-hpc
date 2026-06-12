/**
 *
 * @file openmp/codelet_dzasum.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dzasum OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, int M, int N,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

#pragma omp task firstprivate( storev, uplo, M, N, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    TCORE_dzasum( storev, uplo, M, N, tileA, tileB->mat );

    (void)options;
}
