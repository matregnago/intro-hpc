/**
 *
 * @file starpu/runtime_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU IPIV descriptor routines. These routines are used by laswp operations.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @author Samuel Thibault
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"

/**
 *  Create ws_pivot runtime structures
 */
void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv )
{
    assert( ipiv );
    size_t                nbhandles = 3 * ipiv->mt;
    starpu_data_handle_t *handles   = calloc( nbhandles, sizeof(starpu_data_handle_t) );
    ipiv->ipiv = handles;
    handles   += ipiv->mt;
    ipiv->perm = handles;
    handles   += ipiv->mt;
    ipiv->invp = handles;
#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe pivot structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        ipiv->mpitag_ipiv = chameleon_starpu_tag_book( nbhandles );
        if ( ipiv->mpitag_ipiv == -1 ) {
            chameleon_fatal_error("RUNTIME_ipiv_create", "Can't pursue computation since no more tags are available for ipiv structure");
            return;
        }
        ipiv->mpitag_perm = ipiv->mpitag_ipiv + ipiv->mt;
        ipiv->mpitag_invp = ipiv->mpitag_perm + ipiv->mt;
    }
#endif
}

/**
 *  Destroy ws_pivot runtime structures
 */
void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    size_t                i;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->ipiv);
    size_t                nbhandles = 3 * ipiv->mt;

    if ( handle ) {
        for(i=0; i<nbhandles; i++) {
            if ( *handle != NULL ) {
                starpu_data_unregister( *handle );
                *handle = NULL;
            }
            handle++;
        }

        free( ipiv->ipiv    );
        ipiv->ipiv    = NULL;
        ipiv->perm    = NULL;
        ipiv->invp    = NULL;
        chameleon_starpu_tag_release( ipiv->mpitag_ipiv );
    }
}

static inline void*
__runtime_ipiv_getaddr( const CHAM_ipiv_t    *ipiv,
                        starpu_data_handle_t *handle,
                        int64_t               tagbase,
                        int                   m,
                        int                  *data )
{
    int     ncols, home_node = -1;
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    if ( tagbase == ipiv->mpitag_ipiv ) {
        ncols = (mm == (ipiv->mt-1)) ? ipiv->m - mm * ipiv->mb : ipiv->mb;
    }
    else {
        if ( ipiv->withidx ) {
            ncols = (ipiv->max_mt - m + 1) + 2 * ipiv->mb;
        }
        else {
            ncols = ipiv->mb;
        }
    }

    if ( data ) {
        data += ipiv->i + m * ipiv->mb;
        home_node = STARPU_MAIN_RAM;
    }
    starpu_vector_data_register( handle, home_node, (uintptr_t)data, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int                owner = ipiv->get_rankof( ipiv, m, m );
        int64_t            tag   = tagbase + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    return __runtime_ipiv_getaddr( ipiv, ipiv->ipiv, ipiv->mpitag_ipiv, m, ipiv->data );
}

void *RUNTIME_ipiv_getperm( const CHAM_ipiv_t *ipiv, int m )
{
    return __runtime_ipiv_getaddr( ipiv, ipiv->perm, ipiv->mpitag_perm, m, NULL );
}

void *RUNTIME_ipiv_getinvp( const CHAM_ipiv_t *ipiv, int m )
{
    return __runtime_ipiv_getaddr( ipiv, ipiv->invp, ipiv->mpitag_invp, m, NULL );
}

static inline void
__runtime_ipiv_flushone( RUNTIME_sequence_t   *sequence,
                         const CHAM_ipiv_t    *ipiv,
                         starpu_data_handle_t *handle,
                         int64_t               shift )
{
    handle += shift;

    if ( *handle == NULL ) {
        return;
    }

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush( sequence->comm, *handle );
    if ( starpu_mpi_data_get_rank( *handle ) == ipiv->myrank )
#endif
    {
        chameleon_starpu_data_wont_use( *handle );
    }

    (void)sequence;
    (void)ipiv;
}

void RUNTIME_ipiv_flushone( RUNTIME_sequence_t *sequence,
                            CHAM_ipiv_e which, const CHAM_ipiv_t *ipiv, int m )
{
    int64_t mm = m + ( ipiv->i / ipiv->mb );

    if ( which & CHAMIPIV_IPIV ) {
        __runtime_ipiv_flushone( sequence, ipiv, ipiv->ipiv, mm );
    }
    if ( which & CHAMIPIV_PERM ) {
        __runtime_ipiv_flushone( sequence, ipiv, ipiv->perm, mm );
    }
    if ( which & CHAMIPIV_INVP ) {
        __runtime_ipiv_flushone( sequence, ipiv, ipiv->invp, mm );
    }
}


void RUNTIME_ipiv_flushall( RUNTIME_sequence_t *sequence,
                            CHAM_ipiv_e which, const CHAM_ipiv_t *ipiv )
{
    int m;

    for (m = 0; m < ipiv->mt; m++)
    {
        RUNTIME_ipiv_flushone( sequence, which, ipiv, m );
    }
}

void RUNTIME_ipiv_gather( RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t  *desc,
                          int                *ipiv,
                          int                 node )
{
    int64_t mt   = desc->mt;
    int64_t mb   = desc->mb;
    int64_t tag  = chameleon_starpu_tag_book( (int64_t)(desc->mt) );
    int     rank = CHAMELEON_Comm_rank();
    int     m;

    for (m = 0; m < mt; m++, ipiv += mb) {
        starpu_data_handle_t ipiv_src = RUNTIME_ipiv_getaddr( desc, m );

#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_get_data_on_node( sequence->comm, ipiv_src, node );
        if ( rank == node )
#endif
        {
            starpu_data_handle_t ipiv_dst;
            int       ncols     = (m == (mt-1)) ? desc->m - m * mb : mb;
            uintptr_t ipivptr   = (rank == node) ? (uintptr_t)ipiv : 0;
            int       home_node = (rank == node) ? STARPU_MAIN_RAM : -1;

            starpu_vector_data_register( &ipiv_dst, home_node, ipivptr, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
            starpu_mpi_data_register( ipiv_dst, tag + m, 0 );
#endif /* defined(CHAMELEON_USE_MPI) */

            assert( ipiv_dst );

            starpu_data_cpy( ipiv_dst, ipiv_src, 0, NULL, NULL );
            starpu_data_unregister( ipiv_dst );
        }
    }

    chameleon_starpu_tag_release( tag );
}

