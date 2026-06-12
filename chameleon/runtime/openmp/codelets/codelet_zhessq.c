/**
 *
 * @file openmp/codelet_zhessq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhessq OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    INSERT_TASK_zsyssq( options, storev, uplo, n,
                        A, Am, An,
                        SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
