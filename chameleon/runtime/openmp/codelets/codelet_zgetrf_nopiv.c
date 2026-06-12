/**
 *
 * @file openmp/codelet_zgetrf_nopiv.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv OpenMP codelet
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

void INSERT_TASK_zgetrf_nopiv( const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    int info = 0;
#pragma omp task firstprivate( m, n, ib, tileA, info ) depend( inout:tileA[0] )
    TCORE_zgetrf_nopiv( m, n, ib, tileA, &info );

    (void)options;
    (void)nb;
    (void)iinfo;
}
