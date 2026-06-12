/**
 *
 * @file openmp/codelet_ztrssq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrssq OpenMP codelet
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

void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_diag_t diag,
                        int m, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *W, int Wm, int Wn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileW = W->get_blktile( W, Wm, Wn );
#pragma omp task firstprivate( uplo, diag, m, n, tileA, tileW ) depend( in:tileA[0] ) depend( inout:tileW[0] )
    TCORE_ztrssq( uplo, diag, m, n, tileA, tileW );

    (void)options;
}
