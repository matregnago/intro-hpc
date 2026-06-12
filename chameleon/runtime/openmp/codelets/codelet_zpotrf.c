/**
 *
 * @file openmp/codelet_zpotrf.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf OpenMP codelet
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
#include "coreblas.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       int iinfo )
{
    int info = 0;
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
#pragma omp task firstprivate( uplo, n, info, tileA ) depend( inout:tileA[0] )
    TCORE_zpotrf( uplo, n, tileA, &info );

    (void)options;
    (void)nb;
    (void)iinfo;
}
