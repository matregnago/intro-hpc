/**
 *
 * @file parsec/runtime_pivot.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC descriptor routines
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
#include "chameleon_parsec.h"

void RUNTIME_pivot_create( CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)pivot;
}

void RUNTIME_pivot_destroy_submit( RUNTIME_sequence_t *sequence,
                                   CHAM_desc_pivot_t  *pivot )
{
    assert( 0 );
    (void)sequence;
    (void)pivot;
}

void RUNTIME_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)pivot;
}

void *RUNTIME_pivot_getaddr( const CHAM_desc_pivot_t *pivot,
                             int rank, int h )
{
    assert( 0 );
    (void)pivot;
    (void)h;
    return NULL;
}

void RUNTIME_pivot_flushone( RUNTIME_sequence_t      *sequence,
                             const CHAM_desc_pivot_t *pivot, int rank )
{
    assert( 0 );
    (void)sequence;
    (void)pivot;
    (void)rank;
}

void RUNTIME_pivot_flushall( RUNTIME_sequence_t      *sequence,
                             const CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)pivot;
    (void)sequence;
}

void RUNTIME_pivot_invalidate( const CHAM_desc_pivot_t *pivot,
                               int                      rank,
                               int                      h )
{
    (void)pivot;
    (void)h;
    (void)rank;
}
