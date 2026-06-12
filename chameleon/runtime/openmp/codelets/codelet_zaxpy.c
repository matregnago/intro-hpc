/**
 *
 * @file openmp/codelet_zaxpy.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy OpenMP codelet
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

void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                      int M, CHAMELEON_Complex64_t alpha,
                      const CHAM_desc_t *A, int Am, int An, int incA,
                      const CHAM_desc_t *B, int Bm, int Bn, int incB )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
#pragma omp task firstprivate( M, alpha, incA, incB, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    TCORE_zaxpy( M, alpha, tileA, incA, tileB, incB );

    (void)options;
}
