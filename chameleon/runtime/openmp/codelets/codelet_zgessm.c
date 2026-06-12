/**
 *
 * @file openmp/codelet_zgessm.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm OpenMP codelet
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

void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       const CHAM_desc_t *L, int Lm, int Ln,
                       const CHAM_desc_t *D, int Dm, int Dn,
                       const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t *tileD = D->get_blktile( D, Dm, Dn );
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

#pragma omp task firstprivate( m, n, k, ib, IPIV, tileD, tileA ) depend( in:tileD[0] ) depend( inout:tileA[0] )
    TCORE_zgessm( m, n, k, ib, IPIV, tileD, tileA );

    (void)options;
    (void)nb;
    (void)L;
    (void)Lm;
    (void)Ln;
}
