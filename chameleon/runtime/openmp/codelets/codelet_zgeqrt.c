/**
 *
 * @file openmp/codelet_zgeqrt.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrt OpenMP codelet
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

void INSERT_TASK_zgeqrt( const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    int ws_size = options->ws_wsize;

#pragma omp task firstprivate( ws_size, m, n, ib, tileA, tileT ) depend( inout:tileA[0] ) depend( out:tileT[0] )
    {
      CHAMELEON_Complex64_t TAU[ws_size];
      CHAMELEON_Complex64_t *work = TAU + chameleon_max( m, n );

      TCORE_zlaset( ChamUpperLower, ib, n, 0., 0., tileT );
      TCORE_zgeqrt( m, n, ib, tileA, tileT, TAU, work );
    }

    (void)options;
    (void)nb;
}
