/**
 *
 * @file openmp/codelet_zsyssq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyssq OpenMP codelet
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

void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                        cham_store_t storev, cham_uplo_t uplo, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileSCALESUMSQ = SCALESUMSQ->get_blktile( SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
#pragma omp task firstprivate( storev, uplo, n, tileA, tileSCALESUMSQ ) depend( in:tileA[0] ) depend( inout:tileSCALESUMSQ[0] )
    TCORE_zsyssq( storev, uplo, n, tileA, tileSCALESUMSQ );

    (void)options;
}
