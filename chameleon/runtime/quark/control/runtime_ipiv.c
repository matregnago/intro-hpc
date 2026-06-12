/**
 *
 * @file quark/runtime_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark descriptor routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Matteo Marcos
 * @date 2025-12-11
 *
 */
#include "chameleon_quark.h"

void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)ipiv;
}

void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)ipiv;
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)ipiv;
    (void)m;
    return NULL;
}

void *RUNTIME_ipiv_getperm( const CHAM_ipiv_t *ipiv, int k )
{
    assert( 0 );
    (void)ipiv;
    (void)k;
    return NULL;
}

void *RUNTIME_ipiv_getinvp( const CHAM_ipiv_t *ipiv, int k )
{
    assert( 0 );
    (void)ipiv;
    (void)k;
    return NULL;
}

void RUNTIME_ipiv_flushone( RUNTIME_sequence_t *sequence,
                            CHAM_ipiv_e         which,
                            const CHAM_ipiv_t  *ipiv,
                            int                 m )
{
    assert( 0 );
    (void)sequence;
    (void)which;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_flushall( RUNTIME_sequence_t *sequence,
                            CHAM_ipiv_e         which,
                            const CHAM_ipiv_t  *ipiv )
{
    assert( 0 );
    (void)sequence;
    (void)which;
    (void)ipiv;
}

void RUNTIME_ipiv_gather( RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t  *desc,
                          int                *ipiv,
                          int                 node )
{
    assert( 0 );
    (void)sequence;
    (void)desc;
    (void)ipiv;
    (void)node;
}
