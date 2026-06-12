/**
 *
 * @file openmp/codelet_zgessq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessq OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                        cham_store_t storev, int m, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileScaleSum = SCALESUMSQ->get_blktile( SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
#pragma omp task firstprivate( storev, m, n, tileA, tileScaleSum ) depend( in:tileA[0] ) depend( inout:tileScaleSum[0] )
    TCORE_zgessq( storev, m, n, tileA, tileScaleSum );

    (void)options;
}
