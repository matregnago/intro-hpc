/**
 *
 * @file openmp/codelet_zlaset2.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 OpenMP codelet
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

void INSERT_TASK_zlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An )
{

    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
#pragma omp task firstprivate( uplo, m, n, alpha, tileA ) depend( inout:tileA[0] )
    TCORE_zlaset2( uplo, m, n, alpha, tileA );

    (void)options;
}
