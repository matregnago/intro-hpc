/**
 *
 * @file starpu/runtime_pivot.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU pivot descriptor routines. These routines are used
 * exclusively by the LU partial pivoting panel kernel.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-11
 *
 */
#include "chameleon_starpu_internal.h"

/**
 *  Create ws_pivot runtime structures
 */
void RUNTIME_pivot_create( CHAM_desc_pivot_t *pivot )
{
    assert( pivot );
    size_t                nbhandles = 2 * pivot->P;
    starpu_data_handle_t *handles   = calloc( nbhandles, sizeof(starpu_data_handle_t) );
    pivot->nextpiv = handles;
    handles += pivot->P;
    pivot->prevpiv = handles;
#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe pivot structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        pivot->mpitag_nextpiv = chameleon_starpu_tag_book( nbhandles );
        if ( pivot->mpitag_nextpiv == -1 ) {
            chameleon_fatal_error("RUNTIME_pivot_create", "Can't pursue computation since no more tags are available for pivot structure");
            return;
        }
        pivot->mpitag_prevpiv = pivot->mpitag_nextpiv + pivot->P;
    }
#endif
}

/**
 *  Asynchronously destroy ws_pivot runtime structures
 */
void RUNTIME_pivot_destroy_submit( RUNTIME_sequence_t *sequence,
                                   CHAM_desc_pivot_t  *pivot )
{
    size_t                i;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(pivot->nextpiv);
    size_t                nbhandles = 2 * pivot->P;

    if ( !handle ) {
        return;
    }

    for ( i = 0; i < nbhandles; i++ ) {
        if ( *handle != NULL ) {
            starpu_data_unregister_submit( *handle );
            *handle = NULL;
        }
        handle++;
    }

    free( pivot->nextpiv );
    pivot->nextpiv = NULL;
    pivot->prevpiv = NULL;
    (void)sequence;

    /*
     * WARNING: tags are not released in submit as they may always been used by
     * the runtime it can only be done by the synchronous destroy
     */
}

/**
 *  Destroy ws_pivot runtime structures
 */
void RUNTIME_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(pivot->nextpiv);

    if ( handle ) {
        size_t i;
        size_t nbhandles = 2 * pivot->P;

        for ( i = 0; i < nbhandles; i++ ) {
            if ( *handle != NULL ) {
                starpu_data_unregister( *handle );
                *handle = NULL;
            }
            handle++;
        }

        free( pivot->nextpiv );
        pivot->nextpiv = NULL;
        pivot->prevpiv = NULL;
    }

    /*
     * WARNING: tags are released here and not in the submit as they may been
     * used after the unregistration submission. It can only be done by the
     * synchronous destroy.
     */
    chameleon_starpu_tag_release( pivot->mpitag_nextpiv );
}

void *
RUNTIME_pivot_getaddr( const CHAM_desc_pivot_t *pivot, int rank, int h )
{
    starpu_data_handle_t *handle;
    int64_t               tag;
    int                   Q = pivot->Q;

    if ( h%2 == 0 ) {
        handle = (starpu_data_handle_t*)(pivot->nextpiv);
        tag    = pivot->mpitag_nextpiv;
    }
    else {
        handle = (starpu_data_handle_t*)(pivot->prevpiv);
        tag    = pivot->mpitag_prevpiv;
    }

    assert( handle );
    handle += rank / Q;

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    tag += rank / Q;
    cppi_register( handle, pivot->dtyp, pivot->nb, tag, rank );

    assert( *handle );
    return (void*)(*handle);
}

void RUNTIME_pivot_flushone( RUNTIME_sequence_t      *sequence,
                             const CHAM_desc_pivot_t *pivot, int rank )
{
    starpu_data_handle_t *handle;
    int                   Q = pivot->Q;

    handle = (starpu_data_handle_t*)(pivot->nextpiv);
    handle += rank/Q;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    handle = (starpu_data_handle_t*)(pivot->prevpiv);
    handle += rank/Q;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)pivot;
    (void)rank;
}

void RUNTIME_pivot_flushall( RUNTIME_sequence_t      *sequence,
                             const CHAM_desc_pivot_t *pivot )
{
    int m;

    for (m = 0; m < pivot->Q; m++)
    {
        RUNTIME_pivot_flushone( sequence, pivot, m );
    }
}

void RUNTIME_pivot_invalidate( const CHAM_desc_pivot_t *pivot,
                               int                      rank,
                               int                      h )
{
    /* Protection against incorrect h values */
    if ( h < 0 ) {
        return;
    }
    starpu_data_invalidate_submit( RUNTIME_pivot_getaddr( pivot, rank, h ) );
}
