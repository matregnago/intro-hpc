/**
 *
 * @file openmp/codelet_zlascal.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal OpenMP codelet
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

void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha,
                         const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    if ( alpha == 0. ) {
#pragma omp task firstprivate( uplo, m, n, alpha, tileA ) depend( inout:tileA[0] )
        TCORE_zlaset( uplo, m, n, alpha, alpha, tileA );
    }
    else if ( alpha == 1. ) {
        return;
    }
    else {
#pragma omp task firstprivate( uplo, m, n, alpha, tileA ) depend( inout:tileA[0] )
        TCORE_zlascal( uplo, m, n, alpha, tileA );
    }

    (void)options;
    (void)nb;
}
