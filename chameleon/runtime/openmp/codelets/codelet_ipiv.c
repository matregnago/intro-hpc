/**
 *
 * @file openmp/codelet_ipiv.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP codelets to convert pivot to permutations
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks.h"
#include "coreblas.h"

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_desc_pivot_t *pivot, int k, int h, int rank )
{
    assert( 0 );
    (void)options;
    (void)pivot;
    (void)k;
    (void)h;
    (void)rank;
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2, int mt,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    int *ipiv = NULL; // get pointer from ipivdesc
    int *perm = NULL; // get pointer from ipivdesc
    int *invp = NULL; // get pointer from ipivdesc

    if ( ipivdesc->withidx ) {
#pragma omp task firstprivate( m0, m, k, K1, K2, mt ) depend( in:ipiv[0] ) depend( inout:perm[0] ) depend( inout:invp[0] )
        {
            int permtmp[ m ];
            int invptmp[ m ];
            CORE_ipiv_to_perm( m0, m, k, K1, K2,
                               ipiv, permtmp, invptmp );
            CORE_perm_to_idx( m0, m, mt,
                              permtmp, invptmp, perm, invp );
        }
    }
    else {
#pragma omp task firstprivate( m0, m, k, K1, K2 ) depend( in:ipiv[0] ) depend( inout:perm[0] ) depend( inout:invp[0] )
        {
            CORE_ipiv_to_perm( m0, m, k, K1, K2, ipiv, perm, invp );
        }
    }
    (void)options;
    (void)ipivk;
}
