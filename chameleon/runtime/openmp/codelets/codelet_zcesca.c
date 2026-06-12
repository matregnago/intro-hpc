/**
 *
 * @file openmp/codelet_zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca OpenMP codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> s d c z
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zcesca( const RUNTIME_option_t *options,
                         int center, int scale, cham_store_t axis,
                         int m, int n, int mt, int nt,
                         const CHAM_desc_t *Gi, int Gim, int Gin,
                         const CHAM_desc_t *Gj, int Gjm, int Gjn,
                         const CHAM_desc_t *G,  int Gm,  int Gn,
                         const CHAM_desc_t *Di, int Dim, int Din,
                         const CHAM_desc_t *Dj, int Djm, int Djn,
                         CHAM_desc_t *A,  int Am,  int An )
{
    CHAM_tile_t *tileGi = Gi->get_blktile( Gi, Gim, Gin );
    CHAM_tile_t *tileGj = Gj->get_blktile( Gj, Gjm, Gjn );
    CHAM_tile_t *tileG  = G->get_blktile( G, Gm, Gn );
    CHAM_tile_t *tileDi = Di->get_blktile( Di, Dim, Din );
    CHAM_tile_t *tileDj = Dj->get_blktile( Dj, Djm, Djn );
    CHAM_tile_t *tileA  = A->get_blktile( A, Am, An );

#pragma omp task firstprivate( center, scale, axis, m, n, mt, nt, tileGi, tileGj, tileG, tileDi, tileDj, tileA ) depend( in:tileGi[0], tileGj[0], tileG[0], tileDi[0], tileDj[0] ) depend( inout:tileA[0] )
    TCORE_zcesca( center, scale, axis, m, n, mt, nt,
                  tileGi, tileGj, tileG, tileDi, tileDj, tileA );

    (void)options;
}
