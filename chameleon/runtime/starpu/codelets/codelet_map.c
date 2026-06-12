/**
 *
 * @file starpu/codelet_map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"

struct cl_map_args_s {
    cham_uplo_t          uplo;
    int                  m, n;
    cham_map_operator_t *op_fcts;
    void                *op_args;
    const CHAM_desc_t   *desc[1];
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)

static inline void
cl_map_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_map_args_s *clargs  = (struct cl_map_args_s *)(t->cl_arg);
    rectask_args_t       *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t     request = RUNTIME_REQUEST_INITIALIZER;
    int                   i;
    int                   ndata   = t->nbuffers;
    cham_map_data_t       data[ndata];

    starpu_cham_rectask_initrequest( t, &request );

    for ( i = 0; i < ndata; i++ ) {
        const CHAM_desc_t *desc = clargs->desc[i];
        CHAM_tile_t       *tile = desc->get_blktile( desc, clargs->m, clargs->n );
        enum starpu_data_access_mode mode = STARPU_TASK_GET_MODE( t, i );

        data[i].access = starpu_to_cham_access( mode );
        data[i].desc   = tile->mat;
    }

    chameleon_pmap( clargs->uplo, t->nbuffers, data, clargs->op_fcts, clargs->op_args,
                    rtargs->sequence, &request );
}

#endif /* CHAMELEON_USE_RECURSIVE_TASKS */

/*
 * Map with a single tile as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_one_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                              clargs->desc[0], tileA );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_one_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                               clargs->desc[0], tileA );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_one_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                              clargs->desc[0], tileA );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_one, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( map_one, cl_map_one_cpu_func, cl_map_one_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_one, cl_map_one_cpu_func, cl_map_one_cuda_func, STARPU_CUDA_ASYNC )
#endif

/*
 * Map with two tiles as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_two_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                              clargs->desc[0], tileA, clargs->desc[1], tileB );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_two_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                               clargs->desc[0], tileA, clargs->desc[1], tileB );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_two_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                              clargs->desc[0], tileA, clargs->desc[1], tileB );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_two, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( map_two, cl_map_two_cpu_func, cl_map_two_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_two, cl_map_two_cpu_func, cl_map_two_cuda_func, STARPU_CUDA_ASYNC )
#endif

/*
 * Map with three tiles as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_three_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                              clargs->desc[0], tileA, clargs->desc[1], tileB,
                              clargs->desc[2], tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_three_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                               clargs->desc[0], tileA, clargs->desc[1], tileB,
                               clargs->desc[2], tileC );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_three_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                              clargs->desc[0], tileA, clargs->desc[1], tileB,
                              clargs->desc[2], tileC );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_three, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
    CODELETS_GPU( map_three, cl_map_three_cpu_func, cl_map_three_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_three, cl_map_three_cpu_func, cl_map_three_cuda_func, STARPU_CUDA_ASYNC )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    struct cl_map_args_s *clargs  = NULL;
    const char           *cl_name = (op_fcts->name == NULL) ? "map" : op_fcts->name;
    int                   exec    = 0;
    int                   i, readonly = 1;
    size_t                clargs_size = 0;
    uint32_t              where       = 0;
    int                   is_rectask  = 1;
    rectask_args_t       *rtargs      = NULL;
    CHAM_tile_t          *tiles[ndata];
    (void)rtargs;

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    for( i=0; i<ndata; i++ ) {
        if ( data[i].access == ChamRW ) {
            CHAMELEON_ACCESS_RW( data[i].desc, m, n );
            readonly = 0;
        }
        else if ( data[i].access == ChamW ) {
            CHAMELEON_ACCESS_W( data[i].desc, m, n );
            readonly = 0;
        }
        else {
            CHAMELEON_ACCESS_R( data[i].desc, m, n );
        }
    }
    exec = __chameleon_need_exec;
    /* Force execution for read-only functions */
    if ( readonly && __chameleon_need_submit ) {
        exec = 1;
    }
    CHAMELEON_END_ACCESS_DECLARATION;

    for( i=0; i<ndata; i++ ) {
        tiles[i] = (data[i].desc)->get_blktile( data[i].desc, m, n );
 #if defined(CHAMELEON_USE_RECURSIVE_TASKS)
        if ( !(tiles[i]->format & CHAMELEON_TILE_DESC) ) {
            is_rectask = 0;
        }
#else
        is_rectask = 0;
#endif
    }

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + sizeof(CHAM_tile_t*) * (ndata-1) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority;
        memcpy( rtargs->tiles, tiles, ndata * sizeof(CHAM_tile_t*) );
    }
#endif

    if ( is_rectask || exec ) {
        clargs_size = sizeof( struct cl_map_args_s ) + sizeof( CHAM_desc_t * ) * (ndata - 1);
        clargs = malloc( clargs_size );
        clargs->uplo    = uplo;
        clargs->m       = m;
        clargs->n       = n;
        clargs->op_fcts = op_fcts;
        clargs->op_args = op_args;
        for( i=0; i<ndata; i++ ) {
            clargs->desc[i] = data[i].desc;
        }
    }

    /* Refine name */
    for( i=0; i<ndata; i++ ) {
        cl_name = chameleon_codelet_name( cl_name, 1, tiles[i] );
    }

    /* Where to execute */
    if ( op_fcts->cpufunc ) {
        where |= STARPU_CPU;
    }
    if ( op_fcts->cudafunc ) {
        where |= STARPU_CUDA;
    }
    if ( op_fcts->hipfunc ) {
        where |= STARPU_HIP;
    }

    /* Insert the task */
    switch( ndata ) {
    case 1:
        /* Insert the task */
        rt_starpu_insert_task(
            &cl_map_one,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),

            /* Common task arguments */
            INSERT_TASK_COMMON_TASK_PARAMS( map_one ),
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,

            /* Recursive task management */
            INSERT_TASK_RECTASK_PARAMS( map )
            0 );
        break;

    case 2:
        /* Insert the task */
        rt_starpu_insert_task(
            &cl_map_two,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[1].access ), RTBLKADDR( data[1].desc, ChamByte, m, n ),

            /* Common task arguments */
            INSERT_TASK_COMMON_TASK_PARAMS( map_two ),
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,

            /* Recursive task management */
            INSERT_TASK_RECTASK_PARAMS( map )
            0 );
        break;

    case 3:
        /* Insert the task */
        rt_starpu_insert_task(
            &cl_map_three,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[1].access ), RTBLKADDR( data[1].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[2].access ), RTBLKADDR( data[2].desc, ChamByte, m, n ),

            /* Common task arguments */
            INSERT_TASK_COMMON_TASK_PARAMS( map_three ),
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,

            /* Recursive task management */
            INSERT_TASK_RECTASK_PARAMS( map )
            0 );
        break;
    }
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( map, map_one, map, ndata );
    int i;

    /* Update name if provided */
    cl_name = (op_fcts->name == NULL) ? "map" : op_fcts->name;

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    /*
     * Register the data handles and initialize exchanges if needed
     * The location is based on the first tile location, with the assumption
     * that all X_i(m,n), with X in all descriptors involved, are stored on the
     * same node.
     */
    starpu_cham_exchange_init_params( options, &params,
                                      data[0].desc->get_rankof( data[0].desc, m, n ) );
    for( i=0; i<ndata; i++ ) {
        starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs,
                                                    data[i].desc, m, n,
                                                    cham_to_starpu_access( data[i].access ) );
    }

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        return;
    }

    if ( params.do_execute )
    {
        callback_fct_t      callback = NULL;
        size_t              clargs_size;
        int                 ret;
        struct starpu_task *task = starpu_task_create();
        task->cl    = cl;
        task->where = 0;

        /* Where to execute */
        if ( op_fcts->cpufunc ) {
            task->where |= STARPU_CPU;
        }
        if ( op_fcts->cudafunc ) {
            task->where |= STARPU_CUDA;
        }
        if ( op_fcts->hipfunc ) {
            task->where |= STARPU_HIP;
        }

        /* Set codelet parameters */
        clargs_size = sizeof( struct cl_map_args_s ) + sizeof( CHAM_desc_t * ) * (ndata - 1);
        clargs = malloc( clargs_size );
        clargs->uplo    = uplo;
        clargs->m       = m;
        clargs->n       = n;
        clargs->op_fcts = op_fcts;
        clargs->op_args = op_args;
        for( i=0; i<ndata; i++ ) {
            clargs->desc[i] = data[i].desc;
        }

        task->cl_arg      = clargs;
        task->cl_arg_size = clargs_size;
        task->cl_arg_free = 1;

        switch( ndata ) {
        case 3:
            task->cl = &cl_map_three;
            callback = cl_map_three_callback;
            break;
        case 2:
            task->cl = &cl_map_two;
            callback = cl_map_two_callback;
            break;
        case 1:
            task->cl = &cl_map_one;
            callback = cl_map_one_callback;
        }

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, callback );
        task->synchronous = chameleon_max( task->synchronous, op_fcts->synchronous );

        /* Flops */
        //task->flops = ...;

        /* Refine name */
        for( i=0; i<ndata; i++ ) {
            cl_name = chameleon_codelet_name( cl_name, 1,
                                              (data[i].desc)->get_blktile( data[i].desc, m, n ) );
        }

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_map", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
