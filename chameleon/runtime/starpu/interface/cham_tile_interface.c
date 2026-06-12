/**
 *
 * @file starpu/cham_tile_interface.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon tile interface for StarPU
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Gwenole Lucas
 * @author Samuel Thibault
 * @author Abel Calluaud
 * @author Ana Hourcau
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"
#if defined(CHAMELEON_USE_HMATOSS)
#include "coreblas/hmat.h"

static inline void
cti_hmat_destroy( starpu_cham_tile_interface_t *cham_tile_interface )
{
    switch( cham_tile_interface->flttype ) {
#if defined(CHAMELEON_PREC_Z)
    case ChamComplexDouble:
        hmat_zdestroy( cham_tile_interface->tile.mat );
        break;
#endif
#if defined(CHAMELEON_PREC_C)
    case ChamComplexFloat:
        hmat_cdestroy( cham_tile_interface->tile.mat );
        break;
#endif
#if defined(CHAMELEON_PREC_D)
    case ChamRealDouble:
        hmat_ddestroy( cham_tile_interface->tile.mat );
        break;
#endif
#if defined(CHAMELEON_PREC_S)
    case ChamRealFloat:
        hmat_sdestroy( cham_tile_interface->tile.mat );
        break;
#endif
    default:
        STARPU_ASSERT_MSG( 0, "cti_hmat_destroy(): unknown flttype\n" );
    }
    cham_tile_interface->tile.mat = NULL;
}

static inline size_t
cti_get_hmat_required_size( starpu_cham_tile_interface_t *cham_tile_interface )
{
    size_t size = 0;

    if ( (cham_tile_interface->tile.format & CHAMELEON_TILE_HMAT) &&
         (cham_tile_interface->tile.mat != NULL ) )
    {
        switch( cham_tile_interface->flttype ) {
#if defined(CHAMELEON_PREC_Z)
        case ChamComplexDouble:
            size = hmat_zsize( cham_tile_interface->tile.mat );
            break;
#endif
#if defined(CHAMELEON_PREC_C)
        case ChamComplexFloat:
            size = hmat_csize( cham_tile_interface->tile.mat );
            break;
#endif
#if defined(CHAMELEON_PREC_D)
        case ChamRealDouble:
            size = hmat_dsize( cham_tile_interface->tile.mat );
            break;
#endif
#if defined(CHAMELEON_PREC_S)
        case ChamRealFloat:
            size = hmat_ssize( cham_tile_interface->tile.mat );
            break;
#endif
        default:
            STARPU_ASSERT_MSG( 0, "cti_get_hmat_required_size(cham_tile_interface): unknown flttype\n" );
        }
    }

    return size;
}
#else
static inline size_t
cti_get_hmat_required_size( starpu_cham_tile_interface_t *cham_tile_interface  __attribute__((unused)) ) {
    return 0;
}
#endif

int
cti_handle_get_m( starpu_data_handle_t handle )
{
    CHAM_tile_t *tile = cti_handle_get( handle );
    return tile->m;
}

int
cti_handle_get_n( starpu_data_handle_t handle )
{
    CHAM_tile_t *tile = cti_handle_get( handle );
    return tile->n;
}

static void
cti_init( void *data_interface )
{
    starpu_cham_tile_interface_t *cham_tile_interface = data_interface;
    cham_tile_interface->id = STARPU_CHAM_TILE_INTERFACE_ID;
    cham_tile_interface->allocsize = -1;
}

static void
cti_register_data_handle( starpu_data_handle_t  handle,
                          int                   home_node,
                          void                 *data_interface )
{
    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *) data_interface;
    int node;

    for (node = 0; node < STARPU_MAXNODES; node++)
    {
        starpu_cham_tile_interface_t *local_interface = (starpu_cham_tile_interface_t *)
            starpu_data_get_interface_on_node(handle, node);

        memcpy( local_interface, cham_tile_interface,
                sizeof( starpu_cham_tile_interface_t ) );

        if ( node != home_node )
        {
            local_interface->dev_handle = 0;
            local_interface->tile.mat  = NULL;
            local_interface->tile.ld   = -1;
        }
    }
}

static starpu_ssize_t
cti_allocate_data_on_node( void *data_interface, unsigned node )
{
    uintptr_t addr = 0, handle;
    starpu_cham_tile_interface_t *cham_tile_interface =
        (starpu_cham_tile_interface_t *) data_interface;

    uint32_t ld = cham_tile_interface->tile.m;
    starpu_ssize_t allocated_memory;

    allocated_memory = cham_tile_interface->allocsize;
    if ( allocated_memory <= 0 ) {
        return 0;
    }

    handle = starpu_malloc_on_node( node, allocated_memory );

    if ( !handle ) {
        return -ENOMEM;
    }

    if ( starpu_node_get_kind(node) != STARPU_OPENCL_RAM ) {
        addr = handle;
    }

    /* update the data properly */
#if defined(CHAMELEON_KERNELS_TRACE)
    cham_tile_interface->tile.name  = NULL;
#endif
    cham_tile_interface->tile.mat   = (void*)addr;
    cham_tile_interface->tile.ld    = ld;
    cham_tile_interface->dev_handle = handle;

    return allocated_memory;
}

static void
cti_free_data_on_node( void *data_interface, unsigned node )
{
    starpu_cham_tile_interface_t *cham_tile_interface =
        (starpu_cham_tile_interface_t *) data_interface;

#if defined(CHAMELEON_USE_HMATOSS)
    if ( (cham_tile_interface->tile.format & CHAMELEON_TILE_HMAT) &&
         (cham_tile_interface->tile.mat != NULL ) )
    {
        cti_hmat_destroy( cham_tile_interface );
    }
    else
#endif
    {
        assert( (uintptr_t)(cham_tile_interface->tile.mat) == cham_tile_interface->dev_handle );
    }

    starpu_free_on_node( node, cham_tile_interface->dev_handle, cham_tile_interface->allocsize );
    cham_tile_interface->tile.mat = NULL;
    cham_tile_interface->dev_handle = 0;
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static void
cti_reuse_data_on_node( void *dst_data_interface, const void *cached_interface, unsigned node )
{
    (void)node;
    starpu_cham_tile_interface_t *dst_cham_tile =
        (starpu_cham_tile_interface_t *) dst_data_interface;
    starpu_cham_tile_interface_t *cached_cham_tile =
        (starpu_cham_tile_interface_t *) cached_interface;

    /* update the data properly */
#if defined(CHAMELEON_KERNELS_TRACE)
    dst_cham_tile->tile.name  = NULL;
#endif
    dst_cham_tile->tile.mat   = cached_cham_tile->tile.mat;
    dst_cham_tile->tile.ld    = dst_cham_tile->tile.m;
    dst_cham_tile->dev_handle = cached_cham_tile->dev_handle;
}
#endif

static void *
cti_to_pointer( void *data_interface, unsigned node )
{
    (void) node;
    starpu_cham_tile_interface_t *cham_tile_interface = data_interface;

    return (void*)(cham_tile_interface->tile.mat);
}

static size_t
cti_get_size(starpu_data_handle_t handle)
{
    starpu_cham_tile_interface_t *cham_tile_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    size_t elemsize = CHAMELEON_Element_Size( cham_tile_interface->flttype );

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( cham_tile_interface->id == STARPU_CHAM_TILE_INTERFACE_ID,
                       "Error. The given data is not a cham_tile." );
#endif

    return  elemsize * (size_t)(cham_tile_interface->tile.m) * (size_t)(cham_tile_interface->tile.n);
}

static size_t
cti_get_alloc_size(starpu_data_handle_t handle)
{
    starpu_cham_tile_interface_t *cham_tile_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( cham_tile_interface->id == STARPU_CHAM_TILE_INTERFACE_ID,
                       "Error. The given data is not a cham_tile." );
#endif

    STARPU_ASSERT_MSG( cham_tile_interface->allocsize != (size_t)-1,
                       "The cham_tile allocation size needs to be defined" );

    return cham_tile_interface->allocsize;
}

static uint32_t
cti_footprint( starpu_data_handle_t handle )
{
    CHAM_tile_t *tile = cti_handle_get( handle );
    return starpu_hash_crc32c_be( tile->m, tile->n );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static uint32_t
cti_alloc_footprint( starpu_data_handle_t handle )
{
    return starpu_hash_crc32c_be( cti_handle_get_allocsize(handle), 0 );
}
#endif

static int
cti_compare( void *data_interface_a, void *data_interface_b )
{
    starpu_cham_tile_interface_t *cham_tile_a = (starpu_cham_tile_interface_t *) data_interface_a;
    starpu_cham_tile_interface_t *cham_tile_b = (starpu_cham_tile_interface_t *) data_interface_b;

    /* Two matrices are considered compatible if they have the same size */
    return ( cham_tile_a->tile.m  == cham_tile_b->tile.m  )
        && ( cham_tile_a->tile.n  == cham_tile_b->tile.n  )
        && ( cham_tile_a->flttype == cham_tile_b->flttype );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static int
cti_alloc_compare(void *data_interface_a, void *data_interface_b)
{
    starpu_cham_tile_interface_t *cham_tile_a = (starpu_cham_tile_interface_t *) data_interface_a;
    starpu_cham_tile_interface_t *cham_tile_b = (starpu_cham_tile_interface_t *) data_interface_b;

    /* Two matrices are considered compatible if they have the same allocated size */
    return ( cham_tile_a->allocsize == cham_tile_b->allocsize );
}
#endif

static void
cti_display( starpu_data_handle_t handle, FILE *f )
{
    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node(handle, STARPU_MAIN_RAM);

    fprintf( f, "%u\t%u\t",
             cham_tile_interface->tile.m,
             cham_tile_interface->tile.n );
}

static int
cti_pack_data_fullrank( starpu_cham_tile_interface_t *cham_tile_interface,
                        void *ptr )
{
    char *matrix = (void *)cham_tile_interface->tile.mat;

    if ( cham_tile_interface->tile.m == cham_tile_interface->tile.ld ) {
        memcpy( ptr, matrix, cham_tile_interface->allocsize );
    }
    else {
        int   n;
        char *tmpptr = ptr;

        for(n=0; n<cham_tile_interface->tile.n; n++)
        {
            size_t elemsize = CHAMELEON_Element_Size( cham_tile_interface->flttype );
            size_t size = elemsize * (size_t)(cham_tile_interface->tile.m);
            memcpy( tmpptr, matrix, size );
            tmpptr += size;
            matrix += elemsize * (size_t)(cham_tile_interface->tile.ld);
        }
    }
    return 0;
}

static int
cti_pack_data_hmat( starpu_cham_tile_interface_t *cham_tile_interface,
                    void *ptr )
{
#if !defined(CHAMELEON_USE_HMATOSS)
    assert( 0 );
    (void)cham_tile_interface;
    (void)ptr;
#else
    hmat_matrix_t *mat = cham_tile_interface->tile.mat;
    STARPU_ASSERT_MSG( mat != NULL, "cti_pack_data_hmat: Try to pack a NULL pointer\n" );
    switch( cham_tile_interface->flttype ) {
#if defined(CHAMELEON_PREC_Z)
    case ChamComplexDouble:
        hmat_zwrite( mat, ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_C)
    case ChamComplexFloat:
        hmat_cwrite( mat, ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_D)
    case ChamRealDouble:
        hmat_dwrite( mat, ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_S)
    case ChamRealFloat:
        hmat_swrite( mat, ptr );
        break;
#endif
    default:
        STARPU_ASSERT_MSG( 0, "cti_pack_data_hmat: unknown flttype\n" );
    }
#endif
    return 0;
}

static int
cti_pack_data( starpu_data_handle_t handle, unsigned node, void **ptr, starpu_ssize_t *count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node(handle, node);
    size_t size;

    size   = (starpu_ssize_t)(cham_tile_interface->allocsize);
    size  += cti_get_hmat_required_size( cham_tile_interface );
    *count = size + sizeof(size_t) + sizeof(CHAM_tile_t);

    if ( ptr != NULL )
    {
        char *tmp;
        *ptr = (void *)starpu_malloc_on_node_flags( node, *count, 0 );
        tmp = (char*)(*ptr);

        /* Start by the size to allocate on reception */
        memcpy( tmp, &size, sizeof(size_t) );
        tmp += sizeof(size_t);

        /* Copy the tile metadata */
        memcpy( tmp, &(cham_tile_interface->tile), sizeof(CHAM_tile_t) );
        tmp += sizeof(CHAM_tile_t);

        /* Pack the real data */
        if ( cham_tile_interface->tile.format & CHAMELEON_TILE_FULLRANK ) {
            cti_pack_data_fullrank( cham_tile_interface, tmp );
        }
        else if ( cham_tile_interface->tile.format & CHAMELEON_TILE_HMAT ) {
            cti_pack_data_hmat( cham_tile_interface, tmp );
        }
        else {
            STARPU_ASSERT_MSG( 1, "Unsupported format for pack." );
        }
    }

    return 0;
}

static int
cti_unpack_data_fullrank( starpu_cham_tile_interface_t *cham_tile_interface,
                          void *ptr )
{
    char *matrix = (void *)cham_tile_interface->tile.mat;
    assert( cham_tile_interface->tile.format & CHAMELEON_TILE_FULLRANK );
    assert( matrix != NULL );

    if ( cham_tile_interface->tile.m == cham_tile_interface->tile.ld ) {
        memcpy( matrix, ptr, cham_tile_interface->allocsize );
    }
    else {
        int   n;
        char *tmpptr = ptr;

        for(n=0 ; n<cham_tile_interface->tile.n; n++)
        {
            size_t elemsize = CHAMELEON_Element_Size( cham_tile_interface->flttype );
            size_t size = elemsize * (size_t)(cham_tile_interface->tile.m);
            memcpy( matrix, tmpptr, size );
            tmpptr += size;
            matrix += elemsize * (size_t)(cham_tile_interface->tile.ld);
        }
    }
    return 0;
}

static int
cti_unpack_data_hmat( starpu_cham_tile_interface_t *cham_tile_interface,
                      void *ptr )
{
    assert( cham_tile_interface->tile.format & CHAMELEON_TILE_HMAT );
#if !defined(CHAMELEON_USE_HMATOSS)
    assert( 0 );
    (void)cham_tile_interface;
    (void)ptr;
#else
    hmat_matrix_t *mat = NULL;
    switch( cham_tile_interface->flttype ) {
#if defined(CHAMELEON_PREC_Z)
    case ChamComplexDouble:
        mat = hmat_zread( ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_C)
    case ChamComplexFloat:
        mat = hmat_cread( ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_D)
    case ChamRealDouble:
        mat = hmat_dread( ptr );
        break;
#endif
#if defined(CHAMELEON_PREC_S)
    case ChamRealFloat:
        mat = hmat_sread( ptr );
        break;
#endif
    default:
        STARPU_ASSERT_MSG( 0, "cti_unpack_data_hmat: unknown flttype\n" );
    }
    cham_tile_interface->tile.mat = mat;
#endif
    return 0;
}

static int
cti_peek_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node(handle, node);

    char *tmp = ptr;

#if defined(CHAMELEON_USE_MPI_DATATYPES)
    /*
     * We may end up here if an early reception occured before the handle of the
     * received data has been registered. Thus, datatype was not existant and we
     * need to unpack the data ourselves
     */
    STARPU_ASSERT( count == cham_tile_interface->allocsize );
    STARPU_ASSERT( cham_tile_interface->tile.format & CHAMELEON_TILE_FULLRANK );
#else
    {
        CHAM_tile_t dsttile;
        size_t size;

        /* Extract the size of the information to unpack */
        memcpy( &size, tmp, sizeof(size_t) );
        tmp += sizeof(size_t);

        /* Extract the tile metadata of the remote tile */
        memcpy( &dsttile, tmp, sizeof(CHAM_tile_t) );
        tmp += sizeof(CHAM_tile_t);

        assert( ( (dsttile.format & CHAMELEON_TILE_HMAT) && (cham_tile_interface->allocsize == 0   )) ||
                (!(dsttile.format & CHAMELEON_TILE_HMAT) && (cham_tile_interface->allocsize == size)) );

        /*
         * Update with the local information. Data is packed now, and do not
         * need leading dimension anymore
         */
        cham_tile_interface->tile.format = dsttile.format;
        cham_tile_interface->tile.ld = cham_tile_interface->tile.m;

        STARPU_ASSERT( cham_tile_interface->tile.m == dsttile.m );
        STARPU_ASSERT( cham_tile_interface->tile.n == dsttile.n );
        STARPU_ASSERT( count == cham_tile_interface->allocsize + sizeof(size_t) + sizeof(CHAM_tile_t) );
    }
#endif

    /* Unpack the real data */
    if ( cham_tile_interface->tile.format & CHAMELEON_TILE_FULLRANK ) {
        cti_unpack_data_fullrank( cham_tile_interface, tmp );
    }
    else if ( cham_tile_interface->tile.format & CHAMELEON_TILE_HMAT ) {
        cti_unpack_data_hmat( cham_tile_interface, tmp );
    }
    else {
        STARPU_ASSERT_MSG( 1, "Unsupported format for unpack." );
    }

    return 0;
}

static int
cti_unpack_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    cti_peek_data( handle, node, ptr, count );

    /* Free the received information */
    starpu_free_on_node_flags( node, (uintptr_t)ptr, count, 0 );

    return 0;
}

static starpu_ssize_t
cti_describe( void *data_interface, char *buf, size_t size )
{
    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *) data_interface;
#if defined(CHAMELEON_KERNELS_TRACE)
    return snprintf( buf, size, "M%ux%ux%u %s",
                     (unsigned) cham_tile_interface->tile.m,
                     (unsigned) cham_tile_interface->tile.n,
                     (unsigned) cham_tile_interface->flttype,
                     cham_tile_interface->tile.name );
#else
    return snprintf( buf, size, "M%ux%ux%u",
                     (unsigned) cham_tile_interface->tile.m,
                     (unsigned) cham_tile_interface->tile.n,
                     (unsigned) cham_tile_interface->flttype );
#endif
}

static int cti_copy_any_to_any( void *src_interface, unsigned src_node,
                                void *dst_interface, unsigned dst_node, void *async_data )
{
    starpu_cham_tile_interface_t *cham_tile_src = (starpu_cham_tile_interface_t *) src_interface;
    starpu_cham_tile_interface_t *cham_tile_dst = (starpu_cham_tile_interface_t *) dst_interface;
    size_t elemsize = CHAMELEON_Element_Size( cham_tile_src->flttype );
    size_t m        = (size_t)(cham_tile_src->tile.m);
    size_t n        = (size_t)(cham_tile_src->tile.n);
    size_t ld_src   = (size_t)(cham_tile_src->tile.ld);
    size_t ld_dst   = (size_t)(cham_tile_dst->tile.ld);
    int    ret      = 0;

    void *src_mat = CHAM_tile_get_ptr( &(cham_tile_src->tile) );
    void *dst_mat = CHAM_tile_get_ptr( &(cham_tile_dst->tile) );

    assert( ld_src >= m );
    assert( ld_dst >= m );

    assert( m == (size_t)(cham_tile_dst->tile.m) );
    assert( n == (size_t)(cham_tile_dst->tile.n) );

    CHAMELEON_DEBUG( "starpu",
                     "[ANY->ANY] src(%s, type:%s, m=%d, n=%d, ld=%d, ptr:%p) dest(%s, type:%s, m=%d, n=%d, ld=%d, ptr:%p)\n",
                     cham_tile_src->tile.name, CHAM_tile_get_typestr( &(cham_tile_src->tile) ),
                     cham_tile_src->tile.m, cham_tile_src->tile.n, cham_tile_src->tile.ld, src_mat,
                     cham_tile_dst->tile.name, CHAM_tile_get_typestr( &(cham_tile_dst->tile) ),
                     cham_tile_dst->tile.m, cham_tile_dst->tile.n, cham_tile_dst->tile.ld, dst_mat );

    m      = m      * elemsize;
    ld_src = ld_src * elemsize;
    ld_dst = ld_dst * elemsize;
#if defined(HAVE_STARPU_INTERFACE_COPY2D)
    if (starpu_interface_copy2d( (uintptr_t) src_mat, 0, src_node,
                                 (uintptr_t) dst_mat, 0, dst_node,
                                 m, n, ld_src, ld_dst, async_data ) ) {
        ret = -EAGAIN;
    }
#else
    if ( (ld_src == m) && (ld_dst == m) )
    {
        /* Optimize unpartitioned and y-partitioned cases */
        if ( starpu_interface_copy( (uintptr_t) src_mat, 0, src_node,
                                    (uintptr_t) dst_mat, 0, dst_node,
                                    m * n, async_data ) )
        {
            ret = -EAGAIN;
	}
    }
    else
    {
        unsigned y;
        for (y = 0; y < n; y++)
        {
            uint32_t src_offset = y * ld_src;
            uint32_t dst_offset = y * ld_dst;

            if ( starpu_interface_copy( (uintptr_t) src_mat, src_offset, src_node,
                                        (uintptr_t) dst_mat, dst_offset, dst_node,
                                        m, async_data ) )
            {
                ret = -EAGAIN;
            }
        }
    }
#endif

    starpu_interface_data_copy( src_node, dst_node, m * n );

    return ret;
}

static const struct starpu_data_copy_methods cti_copy_methods =
{
    .any_to_any = cti_copy_any_to_any,
};

struct starpu_data_interface_ops starpu_interface_cham_tile_ops =
{
    .init                  = cti_init,
    .register_data_handle  = cti_register_data_handle,
    .allocate_data_on_node = cti_allocate_data_on_node,
    .free_data_on_node     = cti_free_data_on_node,
#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
    .reuse_data_on_node    = cti_reuse_data_on_node,
    .alloc_compare         = cti_alloc_compare,
    .alloc_footprint       = cti_alloc_footprint,
#endif
    .to_pointer            = cti_to_pointer,
    .get_size              = cti_get_size,
    .get_alloc_size        = cti_get_alloc_size,
    .footprint             = cti_footprint,
    .compare               = cti_compare,
    .display               = cti_display,
    .pack_data             = cti_pack_data,
#if defined(HAVE_STARPU_DATA_PEEK)
    .peek_data             = cti_peek_data,
#endif
    .unpack_data           = cti_unpack_data,
    .describe              = cti_describe,
    .copy_methods          =&cti_copy_methods,
    .interfaceid           = STARPU_UNKNOWN_INTERFACE_ID,
    .interface_size        = sizeof(starpu_cham_tile_interface_t),
    .name                  = "STARPU_CHAM_TILE_INTERFACE"
};

void
starpu_cham_tile_register( starpu_data_handle_t *handleptr,
                           int                   home_node,
                           CHAM_tile_t          *tile,
                           cham_flttype_t        flttype )
{
    size_t elemsize = CHAMELEON_Element_Size( flttype );
    starpu_cham_tile_interface_t cham_tile_interface =
        {
            .id         = STARPU_CHAM_TILE_INTERFACE_ID,
            .flttype    = flttype,
            .dev_handle = (intptr_t)(tile->mat),
            .allocsize  = -1,
            .tilesize   = elemsize * (size_t)(tile->m) * (size_t)(tile->n),
        };
    memcpy( &(cham_tile_interface.tile), tile, sizeof( CHAM_tile_t ) );
    /* Overwrite the flttype in case it comes from a data conversion */
    cham_tile_interface.tile.flttype = flttype;

    if ( tile->format & CHAMELEON_TILE_FULLRANK ) {
        cham_tile_interface.allocsize = cham_tile_interface.tilesize;
    }
    else if ( tile->format & CHAMELEON_TILE_DESC ) { /* Needed in case starpu ask for it */
        cham_tile_interface.allocsize = cham_tile_interface.tilesize;
    }
    else if ( tile->format & CHAMELEON_TILE_HMAT ) {
        /* For hmat, allocated data will be handled by hmat library. StarPU cannot allocate it for the library */
        cham_tile_interface.allocsize = 0;
    }

    starpu_data_register( handleptr, home_node, &cham_tile_interface, &starpu_interface_cham_tile_ops );

    switch( cham_get_flttype( flttype ) ) {
#if defined(CHAMELEON_PREC_Z)
    case ChamComplexDouble:
        starpu_data_set_reduction_methods( *handleptr, NULL, &cl_zinit );
        break;
#endif
#if defined(CHAMELEON_PREC_C)
    case ChamComplexFloat:
        starpu_data_set_reduction_methods( *handleptr, NULL, &cl_cinit );
        break;
#endif
#if defined(CHAMELEON_PREC_D)
    case ChamRealDouble:
        starpu_data_set_reduction_methods( *handleptr, NULL, &cl_dinit );
        break;
#endif
#if defined(CHAMELEON_PREC_S)
    case ChamRealFloat:
        starpu_data_set_reduction_methods( *handleptr, NULL, &cl_sinit );
        break;
#endif
    case ChamByte:
        starpu_data_set_reduction_methods( *handleptr, NULL, &cl_binit );
        break;
    default:
        STARPU_ASSERT_MSG( 0, "starpu_cham_tile_register: unsupported flttype\n" );
    }
}

size_t
cti_handle_get_allocsize( starpu_data_handle_t handle )
{
    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( cham_tile_interface->id == STARPU_CHAM_TILE_INTERFACE_ID,
                       "Error. The given data is not a cham_tile." );
#endif

    return cham_tile_interface->allocsize;
}

#if defined(CHAMELEON_USE_MPI_DATATYPES)
int
cti_allocate_datatype_node( starpu_data_handle_t handle,
                            unsigned             node,
                            MPI_Datatype        *datatype )
{
    int ret;

    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node( handle, node );

    size_t m  = cham_tile_interface->tile.m;
    size_t n  = cham_tile_interface->tile.n;
    size_t ld = cham_tile_interface->tile.ld;
    size_t elemsize = CHAMELEON_Element_Size( cham_tile_interface->flttype );

#if defined(CHAMELEON_HAVE_MPI_TYPE_VECTOR_C)
    ret = MPI_Type_vector_c( n, m * elemsize, ld * elemsize, MPI_BYTE, datatype );
#else
    assert( ld * elemsize <= (size_t)INT_MAX );
    ret = MPI_Type_vector( n, m * elemsize, ld * elemsize, MPI_BYTE, datatype );
#endif
    STARPU_ASSERT_MSG(ret == MPI_SUCCESS, "MPI_Type_vector failed");

    ret = MPI_Type_commit( datatype );
    STARPU_ASSERT_MSG(ret == MPI_SUCCESS, "MPI_Type_commit failed");

    return 0;
}

int
cti_allocate_datatype( starpu_data_handle_t handle,
                       MPI_Datatype        *datatype )
{
    return cti_allocate_datatype_node( handle, STARPU_MAIN_RAM, datatype );
}

void
cti_free_datatype( MPI_Datatype *datatype )
{
    MPI_Type_free( datatype );
}
#endif

void
starpu_cham_tile_interface_init()
{
    if ( starpu_interface_cham_tile_ops.interfaceid == STARPU_UNKNOWN_INTERFACE_ID )
    {
        starpu_interface_cham_tile_ops.interfaceid = starpu_data_interface_get_next_id();
#if defined(CHAMELEON_USE_MPI_DATATYPES)
  #if defined(HAVE_STARPU_MPI_INTERFACE_DATATYPE_NODE_REGISTER)
        starpu_mpi_interface_datatype_node_register( starpu_interface_cham_tile_ops.interfaceid,
                                                    cti_allocate_datatype_node,
                                                    cti_free_datatype );
  #else
        starpu_mpi_interface_datatype_register( starpu_interface_cham_tile_ops.interfaceid,
                                                cti_allocate_datatype,
                                                cti_free_datatype );
  #endif
#endif
    }
}

void
starpu_cham_tile_interface_fini()
{
    if ( starpu_interface_cham_tile_ops.interfaceid != STARPU_UNKNOWN_INTERFACE_ID )
    {
#if defined(CHAMELEON_USE_MPI_DATATYPES)
        starpu_mpi_interface_datatype_unregister( starpu_interface_cham_tile_ops.interfaceid );
#endif
    }
}

CHAM_tile_t *
RUNTIME_handle2tile( void *interface )
{
    return cti_interface_get( (starpu_cham_tile_interface_t *)interface );
}
