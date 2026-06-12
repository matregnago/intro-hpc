/**
 *
 * @file openmp/codelet_ztsmlq_hetra1.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq_hetra1 OpenMP codelet
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

void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n,
                              const CHAM_desc_t *A2, int A2m, int A2n,
                              const CHAM_desc_t *V,  int Vm,  int Vn,
                              const CHAM_desc_t *T,  int Tm,  int Tn )
{
    CHAM_tile_t *tileA1 = A1->get_blktile( A1, A1m, A1n );
    CHAM_tile_t *tileA2 = A2->get_blktile( A2, A2m, A2n );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    int ldwork = side == ChamLeft ? ib : nb;
    int ws_size = options->ws_wsize;
#pragma omp task firstprivate( ws_size, side, trans, m1, n1, m2, n2, k, ib, tileA1, tileA2, tileV, tileT, ldwork ) depend( inout:tileA1[0], tileA2[0] ) depend( in:tileT[0], tileV[0] )
    {
      CHAMELEON_Complex64_t work[ws_size];
      TCORE_ztsmlq_hetra1( side, trans, m1, n1, m2, n2, k,
                         ib, tileA1, tileA2, tileV, tileT, work, ldwork );
    }

    (void)options;
    (void)nb;
}
