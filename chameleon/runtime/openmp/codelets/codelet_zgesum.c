/**
 *
 * @file openmp/codelet_zgesum.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesum OpenMP codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zgesum( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SUMS, int SUMSm, int SUMSn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileScaleSum = SUMS->get_blktile( SUMS, SUMSm, SUMSn );
#pragma omp task firstprivate( storev, m, n, tileA, tileScaleSum ) depend( in:tileA[0] ) depend( inout:tileScaleSum[0] )
    TCORE_zgesum( storev, m, n, tileA, tileScaleSum );

    (void)options;
}
