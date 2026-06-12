/**
 *
 * @file starpu/codelet_zperm_allreduce.c
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
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

/**
 * @brief allreduce_init
 *
 */

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
cl_cpui_allreduce_init_cpu_func( void *descr[], void *cl_args )
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

    rows = (CHAMELEON_Complex64_t*)(ws->ws.rows);

    for ( i = 0; i < clargs->k; i++ ) {
        idx   = perm[i] / mb_full;
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

CODELETS_CPU( zperm_allreduce_init, cl_cpui_allreduce_init_cpu_func )

/*
 * Copy the reduce functions into a allreduce structure to be used locally
 */
CODELETS_CPU( zperm_allreduce, cl_cpui_redux_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

/*
 * @brief Initialize ws by copying the lines to be permuted from Wu once all the gets are done.
 */
static inline void
insert_task_zperm_allreduce_init( const RUNTIME_option_t *options,
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

    clargs = malloc( sizeof( struct cl_zperm_init_args_s ) );
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
        &cl_zperm_allreduce_init,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zperm_init_args_s),
        STARPU_R,                   ipiv_handle,
        STARPU_R,                   RTBLKADDR( Wu, ChamComplexDouble, Wum, Wun ),
        STARPU_W,                   RUNTIME_perm_getaddr( ws, Wm, Wn ),
        STARPU_EXECUTE_ON_NODE,     me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
}

/*
 * @brief Task to send ws to another process for the reduction
 */
static inline void
insert_task_zperm_allreduce_send( const RUNTIME_option_t *options,
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
insert_task_zperm_allreduce_recv( const RUNTIME_option_t *options,
                                  CHAM_perm_t            *ws,
                                  int                     me,
                                  int                     src,
                                  int                     m,
                                  int                     n )
{
    void *cpui_handle = ( ws->side == ChamLeft ) ? RUNTIME_perm_getaddr( ws, src, n ) :
                                                   RUNTIME_perm_getaddr( ws, m, src );

    rt_starpu_insert_task(
        &cl_zperm_allreduce,
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
 * @brief Initialize ws by copying the lines to be permuted from Wu once all the gets are done.
 */
static inline void
insert_task_zperm_allreduce_init( const RUNTIME_option_t *options,
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

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_allreduce_init, zperm_allreduce_init, zperm_init, 3 )

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( Wu, ChamComplexDouble, Wum, Wun ), STARPU_R );
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
        chameleon_error( "INSERT_TASK_zperm_allreduce_init", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

/*
 * @brief Task to send ws to another process for the reduction
 */
static inline void
insert_task_zperm_allreduce_send( const RUNTIME_option_t *options,
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
insert_task_zperm_allreduce_recv( const RUNTIME_option_t *options,
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

    INSERT_TASK_COMMON_PARAMETERS( zperm_allreduce, 2 );

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
        chameleon_error( "INSERT_TASK_zperm_allreduce", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, cpui_handle );
    (void)clargs;
    (void)me;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static inline void
zperm_allreduce_chameleon_starpu_task( const RUNTIME_option_t *options,
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
    int  np_iter       = np_involved;
    int  p_recv, p_send, me;
    int  tempkn, tempkm, tempkk;
    int  shift = 1;

    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
        return;
    }

    tempkm = A->get_blkdim( A, ipivk, DIM_m, A->m );
    tempkn = A->get_blkdim( A, ipivk, DIM_n, A->n );
    tempkk = ( ws->side == ChamLeft ) ? tempkm : tempkn;

    insert_task_zperm_allreduce_init( options, dir, A->myrank, tempkk, P, Q, ipiv, ipivk,
                                      Am, An, Wu, Wum, Wun, ws, Wm, Wn );

    /* Get my index in the list */
    for( me = 0; me < np_involved; me++ ) {
        if ( proc_involved[me] == A->myrank ) {
            break;
        }
    }

    /* Submit reduction tree */
    assert( me < np_involved );
    while ( np_iter > 1 ) {
        p_send = proc_involved[ ( me + shift               ) % np_involved ];
        p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ];

        insert_task_zperm_allreduce_send( options, ws, A->myrank, p_send, Wm, Wn );
        insert_task_zperm_allreduce_recv( options, ws, A->myrank, p_recv, Wm, Wn );

        shift   = shift << 1;
        np_iter = chameleon_ceil( np_iter, 2 );
    }
}

void
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
                             cham_dir_t              dir,
                             const CHAM_desc_t      *A,
                             CHAM_desc_t            *U,
                             int                     Um,
                             int                     Un,
                             CHAM_ipiv_t            *ipiv,
                             int                     ipivk,
                             int                     k,
                             int                     n,
                             void                   *ws )
{
    struct chameleon_pzlaswp_s *tmp = (struct chameleon_pzlaswp_s *)ws;
    cham_getrf_allreduce_t      alg = tmp->reduce.alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zperm_allreduce_chameleon_starpu_task( options, dir, A, k, n, ipiv, ipivk,
                                               U, Um, Un, &(tmp->ws), Um, Un, &(tmp->reduce) );
    }
}

void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved )
{
    int p;

    for ( p = 0; p < np; p ++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved )
{
    int   p;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    for ( p = 0; p < np; p++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              ipiv_handle,
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_invp_row( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     k,
                                           int                     n )
{
    int   b, rank;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    for ( b = k+1; (b < A->mt) && ((b-(k+1)) < chameleon_desc_datadist_get_iparam(A, 0)); b ++ ) {
        rank = A->get_rankof( A, b, n );
        if ( rank == A->myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              ipiv_handle,
                                              rank, NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_invp_col( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     m,
                                           int                     k )
{
    int   b, rank;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    for ( b = k+1; (b < A->nt) && ((b-(k+1)) < chameleon_desc_datadist_get_iparam(A, 1)); b ++ ) {
        rank = A->get_rankof( A, m, b );
        if ( rank == A->myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              ipiv_handle,
                                              rank, NULL, NULL );
    }
}

#endif /* defined(CHAMELEON_USE_MPI) */
