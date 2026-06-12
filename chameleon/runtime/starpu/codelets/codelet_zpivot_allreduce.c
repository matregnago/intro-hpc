/**
 *
 * @file starpu/codelet_zpivot_allreduce.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to do the reduction
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)

struct cl_redux_args_s {
    int h;
    int n;
};

#if !defined(CHAMELEON_SIMULATION)
static void
zpivot_allreduce_cpu_func( cppi_interface_t *cppi_me,
                           cppi_interface_t *cppi_src,
                           int               h,
                           int               n )
{
    CHAM_pivot_t          *nextpiv_me  = &(cppi_me->pivot);
    CHAM_pivot_t          *nextpiv_src = &(cppi_src->pivot);
    CHAMELEON_Complex64_t *pivrow_me   = (CHAMELEON_Complex64_t *)(nextpiv_me->pivrow);
    CHAMELEON_Complex64_t *pivrow_src  = (CHAMELEON_Complex64_t *)(nextpiv_src->pivrow);

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout: ");
    cppi_display_dbg( cppi_src, stderr, "Global redux Input: ");

    assert( cppi_me->n         == cppi_src->n         );
    assert( cppi_me->h         == cppi_src->h         );
    assert( cppi_me->flttype   == cppi_src->flttype   );
    assert( cppi_me->arraysize == cppi_src->arraysize );

    if ( cabs( pivrow_src[ h ] ) > cabs( pivrow_me[ h ] ) ) {
        nextpiv_me->blkm0  = nextpiv_src->blkm0;
        nextpiv_me->blkidx = nextpiv_src->blkidx;
        cblas_zcopy( n, pivrow_src, 1, pivrow_me, 1 );
    }

    /* Let's copy the diagonal row if needed */
    if ( ( cppi_src->has_diag == 1 ) &&
         ( cppi_me->has_diag  == -1 ) )
    {
        cblas_zcopy( n, nextpiv_src->diagrow, 1, nextpiv_me->diagrow, 1 );
        /*
         * The interface of the pivot is registered once in starpu so
         * the arraysize is not always correct
         */
        assert( cppi_src->arraysize >= sizeof(CHAMELEON_Complex64_t) * n );
        cppi_me->has_diag = 1;
    }

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout(After): ");
}

static void
cl_zpivot_allreduce_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_args_s *clargs   = (struct cl_redux_args_s *) cl_arg;
    cppi_interface_t       *cppi_me  = ((cppi_interface_t *) descr[0]);
    cppi_interface_t       *cppi_src = ((cppi_interface_t *) descr[1]);
    zpivot_allreduce_cpu_func(  cppi_me, cppi_src, clargs->h, clargs->n );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

CODELETS_CPU( zpivot_allreduce, cl_zpivot_allreduce_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT) /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
INSERT_TASK_zpivot_allreduce_send( const RUNTIME_option_t *options,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     me,
                                   int                     dst,
                                   int                     k,
                                   int                     h )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RUNTIME_pivot_getaddr( pivot, me, h ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    (void)k;
}

static void
INSERT_TASK_zpivot_allreduce_recv( const RUNTIME_option_t *options,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     me,
                                   int                     src,
                                   int                     k,
                                   int                     h,
                                   int                     n )
{
    struct cl_redux_args_s *clargs;
    clargs    = malloc( sizeof( struct cl_redux_args_s ) );
    clargs->h = h;
    clargs->n = n;

    rt_starpu_insert_task(
        &cl_zpivot_allreduce,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_args_s),
        STARPU_RW,                RUNTIME_pivot_getaddr( pivot, me, h ),
        STARPU_R,                 RUNTIME_pivot_getaddr( pivot, src, h ),
        STARPU_EXECUTE_ON_NODE,   me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RUNTIME_pivot_getaddr( pivot, src, h ) );

    (void)k;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
INSERT_TASK_zpivot_allreduce_send( const RUNTIME_option_t *options,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     me,
                                   int                     dst,
                                   int                     k,
                                   int                     h )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zpivot_allreduce_send, 1 )

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_pivot_getaddr( pivot, me, h ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
    (void)k;
}

static void
INSERT_TASK_zpivot_allreduce_recv( const RUNTIME_option_t *options,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     me,
                                   int                     src,
                                   int                     k,
                                   int                     h,
                                   int                     n )
{
    int ret;
    struct starpu_task *task;
    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zpivot_allreduce_recv, zpivot_allreduce, redux, 2 )

    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_pivot_getaddr( pivot, me, h ),
                                                  STARPU_RW );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_pivot_getaddr( pivot, src, h ),
                                                  STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs    = malloc( sizeof( struct cl_redux_args_s ) );
    clargs->h = h;
    clargs->n = n;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_redux_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    /* Refine name */
    task->name = cl_name;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zpivot_allreduce", "Failed to submit the task to StarPU" );
        return;
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, RUNTIME_pivot_getaddr( pivot, src, h ) );

    (void)k;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
zpivot_allreduce_chameleon_starpu_task( const RUNTIME_option_t *options,
                                        CHAM_desc_t            *A,
                                        CHAM_desc_pivot_t      *pivot,
                                        int                    *proc_involved,
                                        int                     k,
                                        int                     h,
                                        int                     n )
{
    int np_involved = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k);
    int np_iter     = np_involved;
    int p_recv, p_send, me;
    int shift = 1;

    if ( h >= pivot->n ) {
        return;
    }

    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
    }
    else {
        for( me = 0; me < np_involved; me++ ) {
            if ( proc_involved[me] == A->myrank ) {
                break;
            }
        }
        assert( me < np_involved );
        while ( np_iter > 1 ) {
            p_send = proc_involved[ ( me + shift               ) % np_involved ];
            p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ];

            INSERT_TASK_zpivot_allreduce_send( options, pivot, A->myrank, p_send, k, h    );
            INSERT_TASK_zpivot_allreduce_recv( options, pivot, A->myrank, p_recv, k, h, n );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}

void
INSERT_TASK_zpivot_allreduce( const RUNTIME_option_t *options,
                             CHAM_desc_t            *A,
                             CHAM_desc_pivot_t      *pivot,
                             int                     k,
                             int                     h,
                             int                     n,
                             void                   *ws )
{
    struct chameleon_pzlaswp_s *tmp = (struct chameleon_pzlaswp_s *)ws;
    cham_getrf_allreduce_t      alg = tmp->reduce.alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zpivot_allreduce_chameleon_starpu_task( options, A, pivot, tmp->reduce.proc_involved, k, h, n );
    }
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
