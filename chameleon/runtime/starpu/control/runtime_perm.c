/**
 *
 * @file starpu/runtime_perm.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU panel permutation update routines. These routines are used by
 * laswp/lapmt operations.
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
#include "chameleon_starpu_internal.h"

void
RUNTIME_perm_create( CHAM_perm_t *ws )
{
    size_t                nbhandles = ( ws->side == ChamLeft ) ? ws->nt * ws->NP :
                                                                 ws->mt * ws->NP;
    starpu_data_handle_t *handles = calloc( nbhandles, sizeof(starpu_data_handle_t) );

    ws->ws = handles;

#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe workspace structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        ws->mpitag_ws = chameleon_starpu_tag_book( nbhandles );
        if ( ws->mpitag_ws == -1 ) {
            chameleon_fatal_error("RUNTIME_perm_create", "Can't pursue computation since no more tags are available for workspace structure");
            return;
        }
    }
#endif
}

void *
RUNTIME_perm_getaddr( const CHAM_perm_t *ws,
                      int                m,
                      int                n )
{
    starpu_data_handle_t *ptr_ws = (starpu_data_handle_t*)(ws->ws);
    int                   ws_idx = ( ws->side == ChamLeft) ? m + n * ws->NP :
                                                             n + m * ws->NP;
    cham_side_t           side   = ws->side;
    int                   ncols, mrows, owner;
    int64_t               tag;

    ptr_ws += ws_idx;
    assert( ptr_ws );

    if ( *ptr_ws != NULL ) {
        return (void*)(*ptr_ws);
    }

    if ( side == ChamLeft ) {
        owner = m;
        ncols = ( n == ws->nt - 1 ) ? ws->n - n * ws->nb : ws->nb ;
        mrows = ws->mb;
    }
    else {
        owner = n;
        ncols = ( m == ws->mt - 1 ) ? ws->m - m * ws->mb : ws->mb ;
        mrows = ws->nb;
    }

    tag = ws->mpitag_ws + ws_idx;

    cpui_register( ptr_ws, side, ws->dtyp, mrows, ncols, tag, owner );

    assert( *ptr_ws );
    return (void*)(*ptr_ws);
}

/**
 *  Destroy workspace runtime structures
 */
void
RUNTIME_perm_destroy( CHAM_perm_t *ws )
{
    size_t                i;
    starpu_data_handle_t *handle    = (starpu_data_handle_t*)(ws->ws);
    size_t                nbhandles = ( ws->side == ChamLeft ) ? ws->nt * ws->NP :
                                                                 ws->mt * ws->NP;

    for( i = 0; i < nbhandles; i++ ) {
        if ( *handle != NULL ) {
            starpu_data_unregister_submit( *handle );
            *handle = NULL;
        }
        handle++;
    }

    free( ws->ws );
    ws->ws = NULL;
    chameleon_starpu_tag_release( ws->mpitag_ws );
}

void
RUNTIME_perm_flush( RUNTIME_sequence_t *sequence,
                    int                 rank,
                    const CHAM_perm_t  *ws,
                    int                 m,
                    int                 n )
{
    starpu_data_handle_t *handle;
    int                   ws_idx = ( ws->side == ChamLeft ) ? m + n * ws->NP :
                                                              n + m * ws->NP;

    handle = (starpu_data_handle_t*)(ws->ws);
    handle += ws_idx;

    if ( *handle == NULL ) {
        return;
    }

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush( sequence->comm, *handle );
    if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
    {
        chameleon_starpu_data_wont_use( *handle );
    }
}
