/**
 *
 * @file openmp/codelet_zssssm.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm OpenMP codelet
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

void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n,
                       const CHAM_desc_t *A2, int A2m, int A2n,
                       const CHAM_desc_t *L1, int L1m, int L1n,
                       const CHAM_desc_t *L2, int L2m, int L2n,
                       const int *IPIV )
{
    CHAM_tile_t *tileA1 = A1->get_blktile( A1, A1m, A1n );
    CHAM_tile_t *tileA2 = A2->get_blktile( A2, A2m, A2n );
    CHAM_tile_t *tileL1 = L1->get_blktile( L1, L1m, L1n );
    CHAM_tile_t *tileL2 = L2->get_blktile( L2, L2m, L2n );

#pragma omp task firstprivate( m1, n1, m2, n2, k, ib, tileA1, tileA2, tileL1, tileL2, IPIV ) \
    depend( inout:tileA1[0], tileA2[0] ) depend( in:tileL1[0], tileL2[0] )
    TCORE_zssssm( m1, n1, m2, n2, k, ib, tileA1, tileA2, tileL1, tileL2, IPIV );

    (void)options;
    (void)nb;
}
