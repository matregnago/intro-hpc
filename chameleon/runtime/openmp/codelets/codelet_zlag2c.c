/**
 *
 * @file openmp/codelet_zlag2c.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_zc.h"
#include "coreblas/coreblas_zctile.h"

void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
#pragma omp task firstprivate( m, n, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        int info = 0;
        TCORE_zlag2c( m, n, tileA, tileB, &info );
    }

    (void)options;
    (void)nb;
}

void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
#pragma omp task firstprivate( m, n, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    TCORE_clag2z( m, n, tileA, tileB );

    (void)options;
    (void)nb;
}
