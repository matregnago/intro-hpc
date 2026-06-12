/**
 *
 * @file starpu/cppi_interface.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon pivot panel interface for StarPU
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"
#undef HAVE_STARPU_REUSE_DATA_ON_NODE

CHAM_pivot_t *
cppi_handle_get( starpu_data_handle_t handle )
{
    cppi_interface_t *cppi_interface = (cppi_interface_t *)
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cppi_interface->id == CPPI_INTERFACE_ID,
                       "Error. The given data is not a CHAM_pivot interface." );
#endif

    return &(cppi_interface->pivot);
}

static void
cppi_init( void *data_interface )
{
    cppi_interface_t *cppi_interface = (cppi_interface_t *)data_interface;
    cppi_interface->id = CPPI_INTERFACE_ID;
    cppi_interface->h  = -1;
    cppi_interface->has_diag = -1;
}

static void
cppi_register_data_handle( starpu_data_handle_t  handle,
                           int                   home_node,
                           void                 *data_interface )
{
    cppi_interface_t *cppi_interface =
        (cppi_interface_t *) data_interface;
    int node;

    for (node = 0; node < STARPU_MAXNODES; node++)
    {
        cppi_interface_t *local_interface = (cppi_interface_t *)
            starpu_data_get_interface_on_node(handle, node);

        memcpy( local_interface, cppi_interface,
                sizeof( cppi_interface_t ) );

        if ( node != home_node )
        {
            local_interface->pivot.pivrow  = NULL;
            local_interface->pivot.diagrow = NULL;
        }
    }
}

static starpu_ssize_t
cppi_allocate_data_on_node( void *data_interface, unsigned node )
{
    cppi_interface_t *cppi_interface = (cppi_interface_t *) data_interface;
    starpu_ssize_t    requested_memory = cppi_interface->arraysize * 2;
    void             *dataptr = NULL;

    dataptr = (void*) starpu_malloc_on_node( node, requested_memory );
    if ( !dataptr ) {
        return -ENOMEM;
    }

    /* WARNING: Should not be a memset if GPU */
    //memset ((void*) dataptr, 0, requested_memory );

    /* update the data properly in consequence */
    cppi_interface->h = -1;
    cppi_interface->has_diag = -1;
    cppi_interface->pivot.pivrow  = dataptr;
    cppi_interface->pivot.diagrow = ((char*)dataptr) + cppi_interface->arraysize;

    return requested_memory;
}

static void
cppi_free_data_on_node( void *data_interface, unsigned node )
{
    cppi_interface_t *cppi_interface   = (cppi_interface_t *) data_interface;
    starpu_ssize_t    requested_memory = cppi_interface->arraysize * 2;

    starpu_free_on_node( node, (uintptr_t)(cppi_interface->pivot.pivrow), requested_memory );

    cppi_interface->pivot.pivrow  = NULL;
    cppi_interface->pivot.diagrow = NULL;
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static void
cppi_reuse_data_on_node( void *dst_data_interface, const void *cached_interface, unsigned node )
{
    (void)node;
    cppi_interface_t *dst_pivot = (cppi_interface_t *)dst_data_interface;
    cppi_interface_t *src_pivot = (cppi_interface_t *)cached_interface;

    /* update the data properly */
    dst_pivot->has_diag = 0;  //src_pivot->has_diag;
    dst_pivot->h        = -1; //src_pivot->h;
    dst_pivot->n        = src_pivot->n;
    dst_pivot->pivot    = src_pivot->pivot;
}
#endif

static size_t
cppi_get_size(starpu_data_handle_t handle)
{
    cppi_interface_t *cppi_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    size_t size;

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cppi_interface->id == CPPI_INTERFACE_ID,
                       "Error. The given data is not a pivot interface." );
#endif

    size = cppi_interface->arraysize * 2 + sizeof(int) * 4;
    return size;
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static size_t
cppi_get_alloc_size(starpu_data_handle_t handle)
{
    cppi_interface_t *cppi_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#if defined(STARPU_DEBUG)
    STARPU_ASSERT_MSG( cppi_interface->id == CPPI_INTERFACE_ID,
                       "Error. The given data is not a pivot interface." );
#endif

    return cppi_interface->arraysize * 2;
}
#endif

static uint32_t
cppi_footprint( starpu_data_handle_t handle )
{
    cppi_interface_t *cppi_interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    return starpu_hash_crc32c_be( 2., cppi_interface->n );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static uint32_t
cppi_alloc_footprint( starpu_data_handle_t handle )
{
    return starpu_hash_crc32c_be( cti_handle_get_allocsize(handle), 0 );
}
#endif

static int
cppi_compare( void *data_interface_a, void *data_interface_b )
{
    cppi_interface_t *cppi_interface_a = (cppi_interface_t *) data_interface_a;
    cppi_interface_t *cppi_interface_b = (cppi_interface_t *) data_interface_b;

    /* Two matrices are considered compatible if they have the same size and the same flttype */
    return ( ( cppi_interface_a->n       == cppi_interface_b->n      ) &&
             ( cppi_interface_a->flttype == cppi_interface_b->flttype) );
}

#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
static int
cppi_alloc_compare(void *data_interface_a, void *data_interface_b)
{
    cppi_interface_t *cppi_a = (cppi_interface_t *) data_interface_a;
    cppi_interface_t *cppi_b = (cppi_interface_t *) data_interface_b;

    /* Two matrices are considered compatible if they have the same allocated size */
    return ( cppi_a->arraysize == cppi_b->arraysize );
}
#endif

static void
cppi_display( starpu_data_handle_t handle, FILE *f )
{
    cppi_interface_t *cppi_interface = (cppi_interface_t *) handle;

    fprintf( f, "%d\t%d\t%d\t%d\n",
             cppi_interface->n,
             cppi_interface->h,
             cppi_interface->pivot.blkm0,
             cppi_interface->pivot.blkidx );
}

static int
cppi_pack_data( starpu_data_handle_t handle, unsigned node, void **ptr, starpu_ssize_t *count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    cppi_interface_t *cppi_interface = (cppi_interface_t *)
        starpu_data_get_interface_on_node(handle, STARPU_MAIN_RAM);

    *count = cppi_get_size( handle );

    if ( ptr != NULL )
    {
        int *tmp;
        *ptr = (void *)starpu_malloc_on_node_flags( node, *count, 0 );
        tmp = (int*)(*ptr);

        /* Copy the tile metadata */
        tmp[0] = cppi_interface->has_diag;
        tmp[1] = cppi_interface->h;
        tmp[2] = cppi_interface->pivot.blkm0;
        tmp[3] = cppi_interface->pivot.blkidx;
        tmp += 4;

        memcpy( tmp, cppi_interface->pivot.pivrow, cppi_interface->arraysize * 2 );
    }

    return 0;
}

static int
cppi_peek_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    STARPU_ASSERT(starpu_data_test_if_allocated_on_node(handle, node));

    cppi_interface_t *cppi_interface = (cppi_interface_t *)
        starpu_data_get_interface_on_node(handle, STARPU_MAIN_RAM);
    int   *tmp  = ptr;
    size_t size = cppi_interface->arraysize * 2;

    STARPU_ASSERT( count == (size + sizeof(int) * 4) );

#if defined(CHAMELEON_USE_MPI_DATATYPES) && 0
    /*
     * We may end up here if an early reception occured before the handle of the
     * received data has been registered. Thus, datatype was not existant and we
     * always need to unpack the data ourselves
     */
#else

    cppi_interface->has_diag     = tmp[0];
    cppi_interface->h            = tmp[1];
    cppi_interface->pivot.blkm0  = tmp[2];
    cppi_interface->pivot.blkidx = tmp[3];
    tmp += 4;

    memcpy( cppi_interface->pivot.pivrow, tmp, size );

#endif
    return 0;
}

static int
cppi_unpack_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    cppi_peek_data( handle, node, ptr, count );

    /* Free the received information */
    starpu_free_on_node_flags( node, (uintptr_t)ptr, count, 0 );

    return 0;
}

static starpu_ssize_t
cppi_describe( void *data_interface, char *buf, size_t size )
{
    cppi_interface_t *cppi_interface = (cppi_interface_t *) data_interface;

    return snprintf( buf, size, "Pivot structure, n %d, h %d, has_diag = %d, blkm0 %d, blkidx %d",
                     cppi_interface->n,
                     cppi_interface->h,
                     cppi_interface->has_diag,
                     cppi_interface->pivot.blkm0,
                     cppi_interface->pivot.blkidx );
}

static int
cppi_copy_any_to_any( void *src_interface, unsigned src_node,
                      void *dst_interface, unsigned dst_node, void *async_data )
{
    cppi_interface_t *cppi_interface_src = (cppi_interface_t *) src_interface;
    cppi_interface_t *cppi_interface_dst = (cppi_interface_t *) dst_interface;
    size_t size;
    int ret = 0;

    STARPU_ASSERT( cppi_interface_src->n       == cppi_interface_dst->n       );
    STARPU_ASSERT( cppi_interface_src->flttype == cppi_interface_dst->flttype );

    cppi_interface_dst->h            = cppi_interface_src->h;
    cppi_interface_dst->has_diag     = cppi_interface_src->has_diag;
    cppi_interface_dst->pivot.blkm0  = cppi_interface_src->pivot.blkm0;
    cppi_interface_dst->pivot.blkidx = cppi_interface_src->pivot.blkidx;

    void *src_mat = cppi_interface_src->pivot.pivrow;
    void *dst_mat = cppi_interface_dst->pivot.pivrow;

    size = cppi_interface_src->arraysize * 2;

    if ( starpu_interface_copy( (uintptr_t) src_mat, 0, src_node,
                                (uintptr_t) dst_mat, 0, dst_node,
                                size, async_data ) )
    {
        ret = -EAGAIN;
    }

    starpu_interface_data_copy( src_node, dst_node, size );

    return ret;
}

static const struct starpu_data_copy_methods cppi_copy_methods =
{
    .any_to_any = cppi_copy_any_to_any,
};

struct starpu_data_interface_ops cppi_ops =
{
    .init                  = cppi_init,
    .register_data_handle  = cppi_register_data_handle,
    .allocate_data_on_node = cppi_allocate_data_on_node,
    .free_data_on_node     = cppi_free_data_on_node,
#if defined(HAVE_STARPU_REUSE_DATA_ON_NODE)
    .reuse_data_on_node    = cppi_reuse_data_on_node,
    .alloc_compare         = cppi_alloc_compare,
    .alloc_footprint       = cppi_alloc_footprint,
#endif
    .get_size              = cppi_get_size,
    .footprint             = cppi_footprint,
    .compare               = cppi_compare,
    .display               = cppi_display,
    .pack_data             = cppi_pack_data,
#if defined(HAVE_STARPU_DATA_PEEK)
    .peek_data             = cppi_peek_data,
#endif
    .unpack_data           = cppi_unpack_data,
    .describe              = cppi_describe,
    .copy_methods          =&cppi_copy_methods,
    .interfaceid           = STARPU_UNKNOWN_INTERFACE_ID,
    .interface_size        = sizeof(cppi_interface_t),
    .name                  = "CPPI_INTERFACE"
};


#if !defined(CHAMELEON_SIMULATION)
static inline int
compare_pivots( cham_flttype_t type, int h, void * redux, void * input ){
    if ( type == ChamRealFloat )
    {
        float *valredux = redux;
        float *valinput = input;
        return fabsf( valredux[h] ) < fabsf( valinput[h] );
    }
    else if ( type == ChamRealDouble )
    {
        double *valredux = redux;
        double *valinput = input;
        return fabs( valredux[h] ) < fabs( valinput[h] );
    }
    else if (type == ChamComplexFloat)
    {
        CHAMELEON_Complex32_t *valredux = redux;
        CHAMELEON_Complex32_t *valinput = input;
        return cabsf( valredux[h] ) < cabsf( valinput[h] );
    }
    else if (type == ChamComplexDouble)
    {
        CHAMELEON_Complex64_t *valredux = redux;
        CHAMELEON_Complex64_t *valinput = input;
        return cabs( valredux[h] ) < cabs( valinput[h] );
    }
    return 0;
}

void
cl_cppi_redux_cpu_func( void *descr[], void *cl_arg )
{
    cppi_interface_t *cppi_redux = ((cppi_interface_t *) descr[0]);
    cppi_interface_t *cppi_input = ((cppi_interface_t *) descr[1]);

    STARPU_ASSERT( cppi_redux->n         == cppi_input->n         );
    STARPU_ASSERT( cppi_redux->flttype   == cppi_input->flttype   );
    STARPU_ASSERT( cppi_redux->arraysize == cppi_input->arraysize );

    cppi_display_dbg( cppi_input, stderr, "BRed Input: ");
    cppi_display_dbg( cppi_redux, stderr, "BRed Inout: ");

    /* Set redux pivot h index to current h index */
    if ( cppi_input->h == -1 ) {
        cppi_input->h = cppi_redux->h;
    }
    if ( cppi_redux->h == -1 ) {
        cppi_redux->h = cppi_input->h;
    }
    assert( cppi_redux->h == cppi_input->h );

    /* Let's copy the diagonal row if needed */
    if ( cppi_input->has_diag == 1 ) {
        assert( cppi_redux->has_diag == -1 );

        memcpy( cppi_redux->pivot.diagrow,
                cppi_input->pivot.diagrow,
                cppi_input->arraysize );
        cppi_redux->has_diag = 1;
    }

    /*
     * Let's now select the pivot:
     * we have to compare the column entry corresponding to the diagonal element.
     */
    {
        int   h            = cppi_redux->h;
        void *pivrow_redux = cppi_redux->pivot.pivrow;
        void *pivrow_input = cppi_input->pivot.pivrow;

        if( compare_pivots( cppi_redux->flttype, h, pivrow_redux, pivrow_input ) )
        {
            cppi_redux->pivot.blkm0  = cppi_input->pivot.blkm0;
            cppi_redux->pivot.blkidx = cppi_input->pivot.blkidx;
            memcpy( pivrow_redux,
                    pivrow_input,
                    cppi_input->arraysize );
        }
    }

    cppi_display_dbg( cppi_input, stderr, "ARed Input: ");
    cppi_display_dbg( cppi_redux, stderr, "ARed Inout: ");

    (void)cl_arg;
    return;
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU(cppi_redux, cl_cppi_redux_cpu_func)

#if !defined(CHAMELEON_SIMULATION)
static void
cl_cppi_init_redux_cpu_func( void *descr[], void *cl_arg )
{
    (void)cl_arg;
    /* (void)descr; */
    cppi_interface_t *cppi_redux = ((cppi_interface_t *) descr[0]);

    /* Redux pivot never has diagonal at initialization */
    cppi_redux->has_diag = -1;
    cppi_redux->h        = -1;

    size_t size = cppi_redux->arraysize;
    memset( cppi_redux->pivot.pivrow,  0, size );
    memset( cppi_redux->pivot.diagrow, 0, size );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( cppi_init_redux, cl_cppi_init_redux_cpu_func )

static void cppi_redux_init( void ) __attribute__( ( constructor ) );
static void cppi_redux_init( void )
{
    cl_cppi_init_redux.nbuffers = 1;
    cl_cppi_init_redux.modes[0] = STARPU_W;
    cl_cppi_init_redux.name = "CPPI ALLREDUX INIT";

    cl_cppi_redux.nbuffers = 2;
    cl_cppi_redux.modes[0] = STARPU_RW | STARPU_COMMUTE;
    cl_cppi_redux.modes[1] = STARPU_R;
    /* Shoulb be RW | COMMUTE to be an allreduce */
    //cl_cppi_redux.modes[1] = STARPU_RW | STARPU_COMMUTE;
    cl_cppi_redux.name = "CPPI ALLREDUX TASK";
}

static void
cppi_set_reduction_methods( starpu_data_handle_t handle)
{
    starpu_data_set_reduction_methods( handle ,
                                       &cl_cppi_redux,
                                       &cl_cppi_init_redux );
}

void
cppi_register( starpu_data_handle_t *handleptr,
               cham_flttype_t        flttype,
               int                   n,
               int64_t               data_tag,
               int                   data_rank )
{
    cppi_interface_t cppi_interface =
        {
            .id = CPPI_INTERFACE_ID,
            .arraysize = n * CHAMELEON_Element_Size( flttype ),
            .flttype = flttype,
            .has_diag = -1,
            .h  = -1,
            .n  = n,
        };
    starpu_data_register( handleptr, -1, &cppi_interface, &cppi_ops );

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_data_register( *handleptr, data_tag, data_rank );
#endif /* defined(CHAMELEON_USE_MPI) */

    cppi_set_reduction_methods( *handleptr );
}

void
cppi_interface_init()
{
    if ( cppi_ops.interfaceid == STARPU_UNKNOWN_INTERFACE_ID )
    {
        cppi_ops.interfaceid = starpu_data_interface_get_next_id();
    }
    cppi_redux_init();
}

void
cppi_interface_fini()
{}
