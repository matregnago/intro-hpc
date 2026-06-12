/**
 *
 * @file starpu/codelet_zperm_reduce_col.c
 *
 * @copyright 2014-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to do the reduction
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)

struct cl_redux_col_args_s {
    int tempnn;
    int nb;
    int m;
    int p;
    int q;
    int me;
};

static void
cl_zperm_reduce_col_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_col_args_s  *clargs     = (struct cl_redux_col_args_s *) cl_arg;
    const CHAM_tile_t           *tileUinout = cti_interface_get( descr[0] );
    const CHAM_tile_t           *tileUin    = cti_interface_get( descr[1] );
    const int                   *perm       = (int *)STARPU_VECTOR_GET_PTR( descr[2] );
    CHAMELEON_Complex64_t       *Uinout     = CHAM_tile_get_ptr( tileUinout );
    const CHAMELEON_Complex64_t *Uin        = CHAM_tile_get_ptr( tileUin );

    int tempnn  = clargs->tempnn;
    int nb      = clargs->nb;
    int m       = clargs->m;
    int p       = clargs->p;
    int q       = clargs->q;
    int me      = clargs->me;
    int mb      = tileUinout->m;
    int i, n, ownerp;

    for ( i = 0; i < tempnn; i++ ) {
        n      = perm[ i ] / nb;
        ownerp = ( (n % q) * p + (m % p) );

        if ( me != ownerp )
        {
            cblas_zcopy( mb, Uin    + i * tileUin->ld,    1,
                             Uinout + i * tileUinout->ld, 1 );
        }
    }
}

CODELETS_CPU( zperm_reduce_col, cl_zperm_reduce_col_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

static void
INSERT_TASK_zperm_reduce_col_send( const RUNTIME_option_t *options,
                                   CHAM_desc_t            *U,
                                   int                     me,
                                   int                     dst,
                                   int                     m )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RTBLKADDR(U, CHAMELEON_Complex64_t, m, me),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
}

static void
INSERT_TASK_zperm_reduce_col_recv( const RUNTIME_option_t *options,
                                   cham_dir_t              dir,
                                   CHAM_desc_t            *U,
                                   CHAM_ipiv_t            *ipiv,
                                   int                     ipivk,
                                   int                     me,
                                   int                     src,
                                   int                     m,
                                   int                     tempnn,
                                   int                     p,
                                   int                     q )
{
    struct cl_redux_col_args_s *clargs;
    void                       *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    clargs = malloc( sizeof( struct cl_redux_col_args_s ) );
    clargs->tempnn  = tempnn;
    clargs->nb      = U->nb;
    clargs->m       = m;
    clargs->p       = p;
    clargs->q       = q;
    clargs->me      = me;

    rt_starpu_insert_task(
        &cl_zperm_reduce_col,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_col_args_s),
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, m, me),
        STARPU_R,                 RTBLKADDR(U, CHAMELEON_Complex64_t, m, src),
        STARPU_R,                 ipiv_handle,
        STARPU_EXECUTE_ON_NODE,   me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, m, src) );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
INSERT_TASK_zperm_reduce_col_send( const RUNTIME_option_t *options,
                                   CHAM_desc_t            *U,
                                   int                     me,
                                   int                     dst,
                                   int                     m )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zperm_reduce_col_send, 1 );

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, m, me ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
}

static void
INSERT_TASK_zperm_reduce_col_recv( const RUNTIME_option_t *options,
                                   cham_dir_t              dir,
                                   CHAM_desc_t            *U,
                                   CHAM_ipiv_t            *ipiv,
                                   int                     ipivk,
                                   int                     me,
                                   int                     src,
                                   int                     m,
                                   int                     tempnn,
                                   int                     p,
                                   int                     q )
{
    int                 ret;
    struct starpu_task *task;
    void               *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_reduce_col_send, zperm_reduce_col, redux_col, 3 );

    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, m, me ),
                                                  STARPU_RW );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, m, src ),
                                                  STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, ipiv_handle, STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_redux_col_args_s ) );
    clargs->tempnn  = tempnn;
    clargs->nb      = U->nb;
    clargs->m       = m;
    clargs->p       = p;
    clargs->q       = q;
    clargs->me      = me;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_redux_col_args_s );
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
        chameleon_error( "INSERT_TASK_zperm_reduce", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, m, src) );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
zperm_reduce_col_chameleon_starpu_task( const RUNTIME_option_t     *options,
                                        cham_dir_t                  dir,
                                        const CHAM_desc_t          *A,
                                        CHAM_desc_t                *U,
                                        int                         Um,
                                        int                         Un,
                                        CHAM_ipiv_t                *ipiv,
                                        int                         ipivk,
                                        int                         m,
                                        int                         k,
                                        CHAM_reduce_t              *reduce )
{
    int *proc_involved = reduce->proc_involved;
    int  np_involved   = chameleon_min( chameleon_desc_datadist_get_iparam(A, 1), A->nt - k );
    int  p_recv, p_send, p/*, me*/;

    if( np_involved == 1 ){
        assert( proc_involved[0] == A->myrank );
        return;
    }

    if ( A->myrank == A->get_rankof( A, m, k ) ) {
        for ( p = 0; p < np_involved; p ++ ) {
            p_recv = proc_involved[p];

            if ( p_recv == A->myrank ) {
                continue;
            }

            INSERT_TASK_zperm_reduce_col_recv( options, dir, U, ipiv, ipivk, A->myrank, p_recv,
                                           m, k == (A->nt-1) ? A->n - k * A->nb : A->nb,
                                           chameleon_desc_datadist_get_iparam(A, 0),
                                           chameleon_desc_datadist_get_iparam(A, 1) );

        }
    }
    else {
        p_send = A->get_rankof( A, m, k );

        INSERT_TASK_zperm_reduce_col_send( options, U, A->myrank, p_send, m );
    }

}

void
INSERT_TASK_zperm_reduce_col( const RUNTIME_option_t *options,
                              cham_dir_t              dir,
                              const CHAM_desc_t      *A,
                              CHAM_desc_t            *U,
                              int                     Um,
                              int                     Un,
                              CHAM_ipiv_t            *ipiv,
                              int                     ipivk,
                              int                     m,
                              int                     k,
                              void                   *ws )
{
    struct chameleon_pzlaswp_s *tmp = (struct chameleon_pzlaswp_s *)ws;
    cham_getrf_allreduce_t      alg = tmp->reduce.alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zperm_reduce_col_chameleon_starpu_task( options, dir, A, U, Um, Un, ipiv, ipivk, m, k, &(tmp->reduce) );
    }
}

#else

void
INSERT_TASK_zperm_reduce_col( const RUNTIME_option_t *options,
                              cham_dir_t              dir,
                              const CHAM_desc_t      *A,
                              CHAM_desc_t            *U,
                              int                     Um,
                              int                     Un,
                              CHAM_ipiv_t            *ipiv,
                              int                     ipivk,
                              int                     m,
                              int                     k,
                              void                   *ws )
{
    (void)options;
    (void)A;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
    (void)ipivk;
    (void)k;
    (void)m;
    (void)ws;
}

#endif
