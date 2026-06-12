/**
 *
 * @file openmp/codelet_zlatro.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro OpenMP codelet
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

void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
#pragma omp task firstprivate( uplo, trans, m, n, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    TCORE_zlatro( uplo, trans, m, n, tileA, tileB );

    (void)options;
    (void)nb;
}
