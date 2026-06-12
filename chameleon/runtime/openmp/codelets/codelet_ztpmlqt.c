/**
 *
 * @file openmp/codelet_ztpmlqt.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmlqt OpenMP codelet
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

void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    int ws_size = options->ws_wsize;

#pragma omp task firstprivate( ws_size, side, trans, M, N, K, L, ib, tileV, tileT, tileA, tileB ) depend( in:tileV[0], tileT[0] ) depend( inout:tileA[0], tileB[0] )
    {
        CHAMELEON_Complex64_t work[ws_size];
        TCORE_ztpmlqt( side, trans, M, N, K, L, ib,
                      tileV, tileT, tileA, tileB, work );
    }

    (void)options;
    (void)nb;
}
