/**
 *
 * @file starpu/cpui_interface.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon panel permutation update interface for StarPU
 *
 * @version 1.4.0
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"
#undef HAVE_STARPU_REUSE_DATA_ON_NODE

static inline size_t
cpui_interface_getoffset( cham_flttype_t flttype, int m )
{
    size_t eltsize = CHAMELEON_Element_Size( flttype );
    return chameleon_ceil( sizeof(int) * m, eltsize ) * eltsize;
}

static size_t
cpui_interface_size( cpui_interface_t *cpui_interface )
{
    int    m       = cpui_interface->m;
    int    n       = cpui_interface->n;
    size_t offset  = cpui_interface->ws.offset;
    size_t eltsize = CHAMELEON_Element_Size( cpui_interface->flttype );
    size_t size    = offset + ( eltsize * n * m );

    return size;
}

CHAM_laswpws_t *
cpui_handle_get( starpu_data_handle_t handle )
{
    cpui_interface_t *cpui_interface = (cpui_interface_t *)
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cpui_interface->id == CPUI_INTERFACE_ID,
                       "Error. The given data is not a CHAM_laswpws interface." );
#endif

    return &(cpui_interface->ws);
}

static void
cpui_init( void *data_interface )
{
    cpui_interface_t *cpui_interface = (cpui_interface_t *)data_interface;
    cpui_interface->id               = CPUI_INTERFACE_ID;
    cpui_interface->n                = -1;
    cpui_interface->m                = -1;
}

static void
cpui_register_data_handle( starpu_data_handle_t  handle,
                           int                   home_node,
                           void                 *data_interface )
{
    cpui_interface_t *cpui_interface =
        (cpui_interface_t *) data_interface;
    int node;

    for ( node = 0; node < STARPU_MAXNODES; node++ ) {
        cpui_interface_t *local_interface = (cpui_interface_t *)
            starpu_data_get_interface_on_node(handle, node);

        memcpy( local_interface, cpui_interface,
                sizeof( cpui_interface_t ) );

        if ( node != home_node ) {
            local_interface->ws.index = NULL;
            local_interface->ws.rows  = NULL;
        }
    }
}

static starpu_ssize_t
cpui_allocate_data_on_node( void     *data_interface,
                            unsigned  node )
{
    cpui_interface_t *cpui_interface = (cpui_interface_t *) data_interface;
    uintptr_t         dataptr        = 0;
    int               offset         = cpui_interface->ws.offset;

    starpu_ssize_t requested_memory_data = cpui_interface_size( cpui_interface );

    dataptr = starpu_malloc_on_node( node, requested_memory_data );
    if ( !dataptr ) {
        return -ENOMEM;
    }

    /* update the data properly in consequence */
    cpui_interface->ws.index = (int*)dataptr;
    cpui_interface->ws.rows  = (void*)(dataptr + offset);

    return requested_memory_data;
}

static void
cpui_free_data_on_node( void     *data_interface,
                        unsigned  node )
{
    cpui_interface_t *cpui_interface        = (cpui_interface_t *) data_interface;
    starpu_ssize_t    requested_memory_data = cpui_interface_size( cpui_interface );

    starpu_free_on_node( node, (uintptr_t)(cpui_interface->ws.index), requested_memory_data );

    cpui_interface->ws.index = NULL;
    cpui_interface->ws.rows  = NULL;
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static void
cpui_reuse_data_on_node( void       *dst_data_interface,
                         const void *cached_interface,
                         unsigned    node )
{
    (void)node;
    cpui_interface_t *dst_ws = (cpui_interface_t *)dst_data_interface;
    cpui_interface_t *src_ws = (cpui_interface_t *)cached_interface;

    /* update the data properly */
    dst_ws->n        = src_ws->n;
    dst_ws->m        = src_ws->m;
    dst_ws->ws.index = src_ws->ws.index;
    dst_ws->ws.rows  = src_ws->ws.rows;
}
#endif

static size_t
cpui_get_size( starpu_data_handle_t handle )
{
    cpui_interface_t *cpui_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    size_t            idx_size, data_size;
    size_t            elt_size = CHAMELEON_Element_Size( cpui_interface->flttype );

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cpui_interface->id == CPUI_INTERFACE_ID,
                       "Error. The given data is not a laswpws interface." );
#endif

    idx_size  = sizeof(int) * cpui_interface->m;
    data_size = elt_size * cpui_interface->ws.nindex * cpui_interface->n;
    return idx_size + data_size;
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static size_t
cpui_get_alloc_size( starpu_data_handle_t handle )
{
    cpui_interface_t *cpui_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cpui_interface->id == CPUI_INTERFACE_ID,
                       "Error. The given data is not a laswpws interface." );
#endif

    return cpui_interface_size( cpui_interface );
}
#endif

static uint32_t
cpui_footprint( starpu_data_handle_t handle )
{
    cpui_interface_t *cpui_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    size_t            size;

    size = cpui_interface_size( cpui_interface );
    return starpu_hash_crc32c_be( 1., size );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static uint32_t
cpui_alloc_footprint( starpu_data_handle_t handle )
{
    return starpu_hash_crc32c_be( cti_handle_get_allocsize(handle), 0 );
}
#endif

static int
cpui_compare( void *data_interface_a,
              void *data_interface_b )
{
    cpui_interface_t *cpui_interface_a = (cpui_interface_t *) data_interface_a;
    cpui_interface_t *cpui_interface_b = (cpui_interface_t *) data_interface_b;

    return ( ( cpui_interface_a->flttype == cpui_interface_b->flttype) &&
             ( cpui_interface_a->m       == cpui_interface_b->m      ) &&
             ( cpui_interface_a->n       == cpui_interface_b->n      ) );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static int
cpui_alloc_compare( void *data_interface_a,
                    void *data_interface_b )
{
    cpui_interface_t *cpui_a = (cpui_interface_t *) data_interface_a;
    cpui_interface_t *cpui_b = (cpui_interface_t *) data_interface_b;

    size_t size_a = cpui_interface_size( cpui_a );
    size_t size_b = cpui_interface_size( cpui_b );

    return ( size_a == size_b );
}
#endif

static void
cpui_display( starpu_data_handle_t  handle,
              FILE                 *f )
{
    cpui_interface_t *cpui_interface = (cpui_interface_t *) handle;

    fprintf( f, "%d\t%d\t%d\t%d\t%d\t",
             cpui_interface->flttype,
             cpui_interface->side,
             cpui_interface->m,
             cpui_interface->n,
             cpui_interface->ws.nindex );
}

static int
cpui_pack_data( starpu_data_handle_t   handle,
                unsigned               node,
                void                 **ptr,
                starpu_ssize_t        *count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    cpui_interface_t *cpui_interface = (cpui_interface_t *)
        starpu_data_get_interface_on_node(handle, STARPU_MAIN_RAM);
    size_t elt_size  = CHAMELEON_Element_Size( cpui_interface->flttype );
    size_t idx_size  = sizeof(int) * cpui_interface->m;
    size_t data_size = elt_size * cpui_interface->ws.nindex * cpui_interface->n;

    *count = sizeof(int) + idx_size + data_size;

    if ( ptr != NULL ) {
        int *tmp;
        *ptr = (void *)starpu_malloc_on_node_flags( node, *count, 0 );
        tmp = (int*)(*ptr);

        /* Copy the tile metadata */
        tmp[0] = cpui_interface->ws.nindex;
        tmp++;

        /* Copy the tile metadata */
        //tmp[0] = cpui_interface->flttype;
        //tmp[1] = cpui_interface->side;
        //tmp[2] = cpui_interface->m;
        //tmp[3] = cpui_interface->n;
        //tmp[4] = cpui_interface->ws.nindex;
        //tmp[5] = cpui_interface->ws.offset;
        //tmp += 6;

        memcpy( tmp, cpui_interface->ws.index, idx_size );
        tmp += cpui_interface->m;

        memcpy( tmp, cpui_interface->ws.rows, data_size );
    }

    return 0;
}

static int
cpui_peek_data( starpu_data_handle_t  handle,
                unsigned              node,
                void                 *ptr,
                size_t                count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    cpui_interface_t *cpui_interface = (cpui_interface_t *)
        starpu_data_get_interface_on_node(handle, STARPU_MAIN_RAM);

    int   *tmp      = ptr;
    size_t elt_size = CHAMELEON_Element_Size( cpui_interface->flttype );
    size_t idx_size = sizeof(int) * cpui_interface->m;
    size_t data_size;

#if defined(CHAMELEON_USE_MPI_DATATYPES) && 0
    /*
     * We may end up here if an early reception occured before the handle of the
     * received data has been registered. Thus, datatype was not existant and we
     * need to unpack the data ourselves
     */
    STARPU_ASSERT( count == (size + sizeof(int) * 4) );

#else

    /* Make sure we at least received the header */
    /* count -= sizeof(int) * 6; */
    /* assert( count >= 0 ); */
    /* cpui_interface->flttype   = tmp[0]; */
    /* cpui_interface->side      = tmp[1]; */
    /* cpui_interface->m         = tmp[2]; */
    /* cpui_interface->n         = tmp[3]; */
    /* cpui_interface->ws.nindex = tmp[4]; */
    /* cpui_interface->ws.offset = tmp[5]; */
    /* tmp += 6; */

    assert( count >= sizeof(int) );
    count -= sizeof(int);
    cpui_interface->ws.nindex = tmp[0];
    tmp++;

    if ( cpui_interface->ws.nindex == 0 ) {
        /* Make sure all indices are set to -1 */
        memset( cpui_interface->ws.index, 0xff, idx_size );
        return 0;
    }

    /* Make sure we received the full index array as annouced */
    assert( count >= idx_size );
    count -= idx_size;

    memcpy( cpui_interface->ws.index, tmp, idx_size );
    tmp += cpui_interface->m;

    data_size = elt_size * cpui_interface->ws.nindex * cpui_interface->n;

    /* Make sure we received the full data array as annouced */
    assert( count >= data_size );
    count -= data_size;

    memcpy( cpui_interface->ws.rows, tmp, data_size );

#endif
    return 0;
}

static int
cpui_unpack_data( starpu_data_handle_t  handle,
                  unsigned              node,
                  void                 *ptr,
                  size_t                count )
{
    cpui_peek_data( handle, node, ptr, count );

    /* Free the received information */
    starpu_free_on_node_flags( node, (uintptr_t)ptr, count, 0 );

    return 0;
}

static starpu_ssize_t
cpui_describe( void   *data_interface,
               char   *buf,
               size_t  size )
{
    cpui_interface_t *cpui_interface = (cpui_interface_t *) data_interface;

    return snprintf( buf, size, "Permutation structure, side %d, n %d, nindex %d",
                     cpui_interface->side,
                     cpui_interface->n,
                     cpui_interface->ws.nindex );
}

static int
cpui_copy_any_to_any( void     *src_interface,
                      unsigned  src_node,
                      void     *dst_interface,
                      unsigned  dst_node,
                      void     *async_data )
{
    cpui_interface_t *cpui_interface_src = (cpui_interface_t *) src_interface;
    cpui_interface_t *cpui_interface_dst = (cpui_interface_t *) dst_interface;
    int               ret = 0;
    void             *src_data,  *dst_data;
    size_t elt_size  = CHAMELEON_Element_Size( cpui_interface_src->flttype );
    size_t idx_size  = sizeof(int) * cpui_interface_src->m;
    size_t data_size = elt_size * cpui_interface_src->ws.nindex * cpui_interface_src->n;

    STARPU_ASSERT( cpui_interface_src->flttype   == cpui_interface_dst->flttype   );
    STARPU_ASSERT( cpui_interface_src->m         == cpui_interface_dst->m         );
    STARPU_ASSERT( cpui_interface_src->n         == cpui_interface_dst->n         );
    STARPU_ASSERT( cpui_interface_src->ws.offset == cpui_interface_dst->ws.offset );

    src_data = cpui_interface_src->ws.index;
    dst_data = cpui_interface_dst->ws.index;

    cpui_interface_dst->side      = cpui_interface_src->side;
    cpui_interface_dst->ws.nindex = cpui_interface_src->ws.nindex;

    data_size = data_size + idx_size;
    if ( starpu_interface_copy( (uintptr_t) src_data, 0, src_node,
                                (uintptr_t) dst_data, 0, dst_node,
                                data_size, async_data ) )
    {
        ret = -EAGAIN;
    }

    starpu_interface_data_copy( src_node, dst_node, data_size );

    return ret;
}

static const struct starpu_data_copy_methods cpui_copy_methods =
{
    .any_to_any = cpui_copy_any_to_any,
};

void
cl_cpui_redux_cpu_func( void *descr[], void *cl_arg )
{
    cpui_interface_t *ws_dst     = (cpui_interface_t*) descr[0];
    cpui_interface_t *ws_src     = (cpui_interface_t*) descr[1];
    char             *rows_src   = (char*)ws_src->ws.rows;
    char             *rows_dst   = (char*)ws_dst->ws.rows;
    int               nindex_dst = ws_dst->ws.nindex;
    int               new_nindex = nindex_dst;
    size_t            eltsize    = CHAMELEON_Element_Size( ws_src->flttype );
    size_t            rowsize    = eltsize * ws_src->n;
    int               ld         = ws_src->n;
    int              *index_src, *index_dst;
    int               src_idx, dst_idx;
    int               i;

    assert( ws_src->n == ws_dst->n );
    assert( ws_src->m == ws_dst->m );

    ld = ws_dst->n;

    index_src = ws_src->ws.index;
    index_dst = ws_dst->ws.index;
    rows_src  = (char*)(ws_src->ws.rows);
    rows_dst  = (char*)(ws_dst->ws.rows);

    for ( i = 0; i < ws_src->m; i++ ){
        /* Skip lines existing in dst buffer */
        if ( index_dst[i] != -1 ) {
            continue;
        }

        /* Skip lines we don't own in src */
        if ( index_src[i] == -1 ) {
            continue;
        }

        src_idx = eltsize * ( index_src[i] * ld );
        dst_idx = eltsize * ( new_nindex   * ld );
        index_dst[ i ] = new_nindex;
        memcpy( rows_dst + dst_idx, rows_src + src_idx, rowsize );
        new_nindex ++;
    }

    ws_dst->ws.nindex = new_nindex;
    assert( new_nindex <= ws_dst->m );

    (void)cl_arg;
}

struct starpu_data_interface_ops cpui_ops =
{
    .init                  = cpui_init,
    .register_data_handle  = cpui_register_data_handle,
    .allocate_data_on_node = cpui_allocate_data_on_node,
    .free_data_on_node     = cpui_free_data_on_node,
#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
    .reuse_data_on_node    = cpui_reuse_data_on_node,
    .alloc_compare         = cpui_alloc_compare,
    .alloc_footprint       = cpui_alloc_footprint,
#endif
    .get_size              = cpui_get_size,
    .footprint             = cpui_footprint,
    .compare               = cpui_compare,
    .display               = cpui_display,
    .pack_data             = cpui_pack_data,
#if defined(HAVE_STARPU_DATA_PEEK)
    .peek_data             = cpui_peek_data,
#endif
    .dontcache             = 1,
    .unpack_data           = cpui_unpack_data,
    .describe              = cpui_describe,
    .copy_methods          =&cpui_copy_methods,
    .interfaceid           = CPUI_INTERFACE_ID,
    .interface_size        = sizeof(cpui_interface_t),
    .name                  = "CPUI_INTERFACE"
};

void
cpui_register( starpu_data_handle_t *handleptr,
               cham_side_t           side,
               cham_flttype_t        flttype,
               int                   m,
               int                   n,
               int64_t               data_tag,
               int                   data_rank )
{
    cpui_interface_t cpui_interface =
        {
            .id      = CPUI_INTERFACE_ID,
            .flttype = flttype,
            .side    = side,
            .m       = m,
            .n       = n,
            .ws      = {
                .index  = NULL,
                .rows   = NULL,
                .nindex = 0,
                .offset = cpui_interface_getoffset( flttype, m ),
            }
        };
    starpu_data_register( handleptr, -1, &cpui_interface, &cpui_ops );

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_data_register( *handleptr, data_tag, data_rank );
#endif /* defined(CHAMELEON_USE_MPI) */
}

void
cpui_interface_init()
{
    if ( cpui_ops.interfaceid == CPUI_INTERFACE_ID )
    {
        cpui_ops.interfaceid = starpu_data_interface_get_next_id();
    }
}

void
cpui_interface_fini()
{}
