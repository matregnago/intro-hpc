/**
 *
 * @file openmp/codelet_zlaset.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

#pragma omp task firstprivate( uplo, m, n, alpha, beta, tileA ) depend( out:tileA[0] )
    TCORE_zlaset( uplo, m, n, alpha, beta, tileA );

    (void)options;
}
