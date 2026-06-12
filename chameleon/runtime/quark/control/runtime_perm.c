/**
 *
 * @file quark/runtime_perm.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark panel permutation update routines. These routines are used by
 * laswp/lapmt operations.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-11
 *
 */
#include "chameleon_quark.h"

void RUNTIME_perm_create( CHAM_perm_t *ws )
{
    assert( 0 );
    (void)ws;
}

void *RUNTIME_perm_getaddr( const CHAM_perm_t *ws,
                            int                m,
                            int                n )
{
    assert( 0 );
    (void)ws;
    (void)m;
    (void)n;
    return NULL;
}

void RUNTIME_perm_destroy( CHAM_perm_t *ws )
{
    assert( 0 );
    (void)ws;
}

void RUNTIME_perm_flush( RUNTIME_sequence_t *sequence,
                          int                 rank,
                          const CHAM_perm_t  *ws,
                          int                 m,
                          int                 n )
{
    assert( 0 );
    (void)sequence;
    (void)ws;
    (void)rank;
    (void)m;
    (void)n;
}

