/**
 *
 * @file openmp/codelet_zgram.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> s d c z
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D,  int Dm,  int Dn,
                              CHAM_desc_t *A,  int Am,  int An )
{
    CHAM_tile_t *tileDi = Di->get_blktile( Di, Dim, Din );
    CHAM_tile_t *tileDj = Dj->get_blktile( Dj, Djm, Djn );
    CHAM_tile_t *tileD  = D->get_blktile( D, Dm, Dn );
    CHAM_tile_t *tileA  = A->get_blktile( A, Am, An );

#pragma omp task firstprivate( uplo, m, n, mt, nt, tileDi, tileDj, tileD, tileA ) depend( in:tileDi[0], tileDj[0], tileD[0] ) depend( inout:tileA[0] )
    TCORE_zgram( uplo, m, n, mt, nt,
                 tileDi, tileDj, tileD, tileA );

    (void)options;
}
