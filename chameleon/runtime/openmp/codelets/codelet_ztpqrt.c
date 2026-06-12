/**
 *
 * @file openmp/codelet_ztpqrt.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    int ws_size = options->ws_wsize;

#pragma omp task firstprivate( ws_size, M, N, L, ib, tileT, tileA, tileB ) depend( inout:tileA[0], tileB[0] ) depend( out:tileT[0] )
    {
      CHAMELEON_Complex64_t tmp[ws_size];

      TCORE_zlaset( ChamUpperLower, ib, N, 0., 0., tileT );
      TCORE_ztpqrt( M, N, L, ib,
                   tileA, tileB, tileT, tmp );
    }

    (void)options;
    (void)nb;
}
