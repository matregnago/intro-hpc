/**
 *
 * @file openmp/codelet_zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP codelets to apply zlaswp on a panel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn )
{
    assert( 0 );
    (void)options;
    (void)side;
    (void)dir;
    (void)m0;
    (void)m;
    (void)n;
    (void)k;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Am;
    (void)An;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t              dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *WA, int WAm, int WAn,
                             const CHAM_desc_t *A,  int Am,  int An )
{
    CHAM_tile_t *tileWA = WA->get_blktile( WA, WAm, WAn );
    CHAM_tile_t *tileA  = A->get_blktile( A, Am, An );
    int         *invp  = NULL; // get invp from ipiv

    assert( tileWA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileA->format  & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( m0, k, ipiv, WA, A ) depend( in:invp ) depend( in:tileWA[0] ) depend( inout:tileA[0] )
    {
        TCORE_zlaswp_set( side, m0, m, n, k, tileWA, tileA, invp );
    }

    (void)options;
    (void)dir;
}

#if defined(CHAMELEON_USE_MPI)

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An )
{
    assert( 0 );
    (void)options;
    (void)ws;
    (void)Wm;
    (void)Wn;
    (void)A;
    (void)Am;
    (void)An;
}

#endif /* if defined(CHAMELEON_USE_MPI) */
