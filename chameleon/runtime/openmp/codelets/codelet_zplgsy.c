/**
 *
 * @file openmp/codelet_zplgsy.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgsy OpenMP codelet
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

void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
#pragma omp task firstprivate( bump, m, n, tileA, bigM, m0, n0, seed ) depend( inout:tileA[0] )
    TCORE_zplgsy( bump, m, n, tileA, bigM, m0, n0, seed );

    (void)options;
}
