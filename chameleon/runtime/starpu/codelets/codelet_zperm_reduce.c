/**
 *
 * @file starpu/codelet_zperm_reduce.c
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
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

struct cl_zperm_init_args_s {
    cham_side_t side;
    int         p;
    int         q;
    int         mb;
    int         Am;
    int         An;
    int         k;
    int         myrank;
};

#if defined(CHAMELEON_USE_MPI)
static void
cl_cpui_reduce_init_cpu_func( void *descr[], void *cl_args )
{

    struct cl_zperm_init_args_s *clargs  = (struct cl_zperm_init_args_s *)cl_args;
    int                         *perm    = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    CHAM_tile_t                 *A_tile  = (CHAM_tile_t *) cti_interface_get( descr[1] );
    cpui_interface_t            *ws      = (cpui_interface_t*) descr[2];
    CHAMELEON_Complex64_t       *A       = CHAM_tile_get_ptr( A_tile );
    int                          ldb     = ws->n;
    cham_side_t                  side    = ws->side;
    int                          mb_full = clargs->mb;
    int                          p       = clargs->p;
    int                          q       = clargs->q;
    int                          myrank  = clargs->myrank;
    int                          Am      = clargs->Am;
    int                          An      = clargs->An;
    CHAMELEON_Complex64_t       *rows;
    int                          nindex;
    int                          i, idx, owner;
    int                          lda, A_inc;

    if ( side == ChamLeft ) {
        A_inc = 1;
        lda   = A_tile->ld;
    }
    else {
        A_inc = A_tile->ld;
        lda   = 1;
    }

    /* Initialize the workspace */
    nindex        = 0;
    ws->ws.nindex = 0;
    memset( ws->ws.index, 0xff, sizeof(int) * ws->m );

    rows = (CHAMELEON_Complex64_t*) ws->ws.rows;

    for ( i = 0; i < clargs->k; i++ ) {
        idx = perm[i] / mb_full;
        owner = ( side == ChamLeft ) ? (idx % p) * q + (An  % q) :
                                       (Am  % p) * q + (idx % q);
        if ( owner != myrank ) {
            continue;
        }

        ws->ws.index[i] = nindex;
        cblas_zcopy( ws->n, A    + i      * A_inc, lda,
                     rows + nindex * ldb,   1 );
        nindex++;
    }
    ws->ws.nindex = nindex;
}

CODELETS_CPU( zperm_reduce_init, cl_cpui_reduce_init_cpu_func )
CODELETS_CPU( zperm_reduce, cl_cpui_redux_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

/*
 * @brief Initialize ws by copying the lines to be permuted from Wu once all the gets are done.
 */
static inline void
insert_task_zperm_reduce_init( const RUNTIME_option_t *options,
                               cham_dir_t              dir,
                               int                     me,
                               int                     k,
                               int                     p,
                               int                     q,
                               CHAM_ipiv_t            *ipiv,
                               int                     ipivk,
                               int                     Am, int An,
                               const CHAM_desc_t      *Wu, int Wum, int Wun,
                               CHAM_perm_t            *ws, int Wm,  int Wn )
{
    struct cl_zperm_init_args_s *clargs;
    void                        *ipiv_handle;

    clargs         = malloc( sizeof( struct cl_zperm_init_args_s ) );
    clargs->side   = ws->side;
    clargs->Am     = Am;
    clargs->An     = An;
    clargs->k      = k;
    clargs->myrank = ipiv->myrank;
    clargs->p      = p;
    clargs->q      = q;
    clargs->mb     = ipiv->mb;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    rt_starpu_insert_task(
        &cl_zperm_reduce_init,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zperm_init_args_s),
        STARPU_R,                 ipiv_handle,
        STARPU_R,                 RTBLKADDR( Wu, ChamComplexDouble, Wum, Wun ),
        STARPU_W,                 RUNTIME_perm_getaddr( ws, Wm, Wn ),
        STARPU_EXECUTE_ON_NODE,   me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
}

/*
 * @brief Task to send ws to another process for the reduction
 */
static inline void
insert_task_zperm_reduce_send( const RUNTIME_option_t *options,
                               CHAM_perm_t            *ws,
                               int                     me,
                               int                     dst,
                               int                     m,
                               int                     n )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RUNTIME_perm_getaddr( ws, m, n ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    (void)me;
}

/*
 * @brief Task to recieve ws from another process and proceed with the reduction.
 */
static inline void
insert_task_zperm_reduce_recv( const RUNTIME_option_t *options,
                               CHAM_perm_t            *ws,
                               int                     me,
                               int                     src,
                               int                     m,
                               int                     n )
{
    void *cpui_handle = ( ws->side == ChamLeft ) ? RUNTIME_perm_getaddr( ws, src, n ) :
                                                   RUNTIME_perm_getaddr( ws, m, src );
    rt_starpu_insert_task(
        &cl_zperm_reduce,
        STARPU_RW,                RUNTIME_perm_getaddr( ws, m, n ),
        STARPU_R,                 cpui_handle,
        STARPU_EXECUTE_ON_NODE,   me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, cpui_handle );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

/*
 * @brief Initialize ws by copying the lines to be permuted from Wu once all the gets are done on Wu.
 */
static inline void
insert_task_zperm_reduce_init( const RUNTIME_option_t *options,
                               cham_dir_t              dir,
                               int                     me,
                               int                     k,
                               int                     p,
                               int                     q,
                               CHAM_ipiv_t            *ipiv,
                               int                     ipivk,
                               int                     Am, int An,
                               const CHAM_desc_t      *Wu, int Wum, int Wun,
                               CHAM_perm_t            *ws, int Wm,  int Wn )
{
    int                 ret;
    struct starpu_task *task;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_reduce_init, zperm_reduce_init, zperm_init, 3 )

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( Wu, ChamComplexDouble, Wum, Wun), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_perm_getaddr( ws, Wm, Wn ), STARPU_W );

    task = starpu_task_create();
    task->cl = cl;

    clargs = malloc( sizeof( struct cl_zperm_init_args_s ) );
    clargs->side   = ws->side;
    clargs->Am     = Am;
    clargs->An     = An;
    clargs->k      = k;
    clargs->myrank = ipiv->myrank;
    clargs->p      = p;
    clargs->q      = q;
    clargs->mb     = ipiv->mb;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zperm_init_args_s );
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
        chameleon_error( "INSERT_TASK_zperm_reduce_init", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

/*
 * @brief Task to send ws to another process for the reduction
 */
static inline void
insert_task_zperm_reduce_send( const RUNTIME_option_t *options,
                               CHAM_perm_t            *ws,
                               int                     me,
                               int                     dst,
                               int                     m,
                               int                     n )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zperm_reduce_send, 1 );

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_perm_getaddr( ws, m, n ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
    (void)me;
}

/*
 * @brief Task to recieve ws from another process and proceed with the reduction.
 */
static inline void
insert_task_zperm_reduce_recv( const RUNTIME_option_t *options,
                               CHAM_perm_t            *ws,
                               int                     me,
                               int                     src,
                               int                     m,
                               int                     n )
{
    int                 ret;
    struct starpu_task *task;
    void               *cpui_handle = ( ws->side == ChamLeft ) ? RUNTIME_perm_getaddr( ws, src, n ) :
                                                                 RUNTIME_perm_getaddr( ws, m, src );

    INSERT_TASK_COMMON_PARAMETERS( zperm_reduce, 2 );

    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_perm_getaddr( ws, m, n ),
                                                  STARPU_RW );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  cpui_handle,
                                                  STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

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
    starpu_mpi_cache_flush( options->sequence->comm, cpui_handle );
    (void)clargs;
    (void)me;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static inline void
zperm_reduce_chameleon_starpu_task( const RUNTIME_option_t *options,
                                    cham_dir_t              dir,
                                    const CHAM_desc_t      *A,
                                    int                     Am,
                                    int                     An,
                                    CHAM_ipiv_t            *ipiv,
                                    int                     ipivk,
                                    const CHAM_desc_t      *Wu,
                                    int                     Wum,
                                    int                     Wun,
                                    CHAM_perm_t            *ws,
                                    int                     Wm,
                                    int                     Wn,
                                    CHAM_reduce_t          *reduce )
{
    int *proc_involved = reduce->proc_involved;
    int  np_involved   = reduce->np_involved;
    int  P             = chameleon_desc_datadist_get_iparam( A, 0 );
    int  Q             = chameleon_desc_datadist_get_iparam( A, 1 );
    int  me            = A->myrank;
    int  arity         = reduce->arity;
    int  root          = A->get_rankof( A, Am, An );
    int  next_np_involved;
    int  step;
    int  me_in_step;
    int  root_in_step;
    int  p_recv, p;
    int  tempkn, tempkm, tempkk;

    if( np_involved == 1 ){
        assert( proc_involved[0] == me );
        return;
    }

    tempkm = A->get_blkdim( A, ipivk, DIM_m, A->m );
    tempkn = A->get_blkdim( A, ipivk, DIM_n, A->n );
    tempkk = ( ws->side == ChamLeft ) ? tempkm : tempkn;

    insert_task_zperm_reduce_init( options, dir, A->myrank, tempkk, P, Q, ipiv, ipivk,
                                   Am, An, Wu, Wum, Wun, ws, Wm, Wn );

    /* Submit reduction tree */
    while ( np_involved != 1 ) {
            next_np_involved = ( np_involved % arity == 0 )
                               ? np_involved / arity
                               : np_involved / arity + 1;

        for ( step = 0; step < next_np_involved; step ++ ) {
            me_in_step = 0;
            root_in_step = 0;

            for ( p = step*arity; ( p < np_involved ) && ( p < (step + 1) * arity ); p ++ ) {
                if ( proc_involved[p] == root ) {
                    root_in_step = 1;
                }
                if ( proc_involved[p] == me ) {
                    me_in_step = 1;
                }
            }

            p_recv = ( root_in_step ) ? root : proc_involved[step * arity];

            if( me == p_recv ) {

                for ( p = step*arity; ( p < np_involved ) && ( p < (step + 1) * arity ); p ++ ) {
                    if ( proc_involved[p] == me ) {
                        continue;
                    }
                    insert_task_zperm_reduce_recv( options, ws, me, proc_involved[p], Wm, Wn);

                }
            }

            else if ( me_in_step ) {

                insert_task_zperm_reduce_send( options, ws, me, p_recv, Wm, Wn );
            }
            proc_involved[step] = p_recv;

        }
        np_involved = next_np_involved;
    }
}

void
INSERT_TASK_zperm_reduce( const RUNTIME_option_t *options,
                          cham_dir_t              dir,
                          const CHAM_desc_t      *A,
                          int                     Am,
                          int                     An,
                          CHAM_ipiv_t            *ipiv,
                          int                     ipivk,
                          const CHAM_desc_t      *Wu,
                          int                     Wum,
                          int                     Wun,
                          void                   *ws,
                          int                     Wm,
                          int                     Wn )
{
    struct chameleon_pzlaswp_s *tmp = (struct chameleon_pzlaswp_s *)ws;
    cham_getrf_allreduce_t      alg = tmp->reduce.alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zperm_reduce_chameleon_starpu_task( options, dir, A, Am, An, ipiv, ipivk, Wu, Wum, Wun, &(tmp->ws), Wm, Wn, &(tmp->reduce) );
    }
}

#else

void
INSERT_TASK_zperm_reduce( const RUNTIME_option_t *options,
                          cham_dir_t              dir,
                          const CHAM_desc_t      *A,
                          int                     m,
                          int                     n,
                          CHAM_ipiv_t            *ipiv,
                          int                     ipivk,
                          const CHAM_desc_t      *Wu,
                          int                     Wum,
                          int                     Wun,
                          void                   *ws,
                          int                     Wm,
                          int                     Wn )
{
    (void)options;
    (void)A;
    (void)m;
    (void)n;
    (void)ws;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)Wu;
    (void)Wum;
    (void)Wun;
    (void)Wm;
    (void)Wn;
}

#endif
