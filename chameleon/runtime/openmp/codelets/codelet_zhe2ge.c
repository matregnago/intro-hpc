/**
 *
 * @file openmp/codelet_zhe2ge.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge OpenMP codelet
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

void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
#pragma omp task firstprivate( uplo, m, n, tileA, tileB ) depend( in: tileA[0] ) depend( inout:tileB[0] )
    TCORE_zhe2ge( uplo, m, n, tileA, tileB );

    (void)options;
    (void)nb;
}
