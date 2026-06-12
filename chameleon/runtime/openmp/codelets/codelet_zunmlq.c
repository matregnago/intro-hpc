/**
 *
 * @file openmp/codelet_zunmlq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq OpenMP codelet
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

void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m, int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );
    int ws_size = options->ws_wsize;
#pragma omp task firstprivate( ws_size, side, trans, m, n, k, ib, nb, tileA, tileT, tileC ) depend( in:tileA[0], tileT[0] ) depend( inout:tileC[0] )
    {
      CHAMELEON_Complex64_t work[ws_size];
      TCORE_zunmlq( side, trans, m, n, k, ib,
                  tileA, tileT, tileC, work, nb );
    }

    (void)options;
    (void)nb;
}
