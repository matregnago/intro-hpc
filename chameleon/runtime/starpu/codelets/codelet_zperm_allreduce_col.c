/**
 *
 * @file starpu/codelet_zperm_allreduce_col.c
 *
 * @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @author Philippe Swartvagher
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
    int p_first;
    int me;
    int shift;
    int np_inv;
};

static void
cl_zperm_allreduce_col_cpu_func( void *descr[], void *cl_arg )
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
    int p_first = clargs->p_first / p;
    int shift   = clargs->shift;
    int np      = clargs->np_inv;
    int me      = ( q <= np ) ? clargs->me / p : ( ( clargs->me / p ) - p_first + q ) % q;
    int mb      = tileUinout->m;
    int first   = me - 2 * shift + 1;
    int last    = me -     shift;
    int i, n, ownerp;

    for ( i = 0; i < tempnn; i++ ) {
        n      = perm[ i ] / nb;
        ownerp = ( q <= np ) ? ( (n % q) * p + (m % p) ) / p : ( ( (n % q) * p + (m % p) ) / p - p_first + q ) % q;

        if ( ( (first    <= ownerp) && (ownerp <= last   ) ) ||
             ( (first+np <= ownerp) && (ownerp <= last+np) ) )
        {
            cblas_zcopy( mb, Uin    + i * tileUin->ld,    1,
                             Uinout + i * tileUinout->ld, 1 );
        }
    }
}

CODELETS_CPU( zperm_allreduce_col, cl_zperm_allreduce_col_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

static void
INSERT_TASK_zperm_allreduce_send_col( const RUNTIME_option_t *options,
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
INSERT_TASK_zperm_allreduce_recv_col( const RUNTIME_option_t *options,
                                      cham_dir_t              dir,
                                      CHAM_desc_t            *U,
                                      CHAM_ipiv_t            *ipiv,
                                      int                     ipivk,
                                      int                     me,
                                      int                     src,
                                      int                     m,
                                      int                     tempnn,
                                      int                     p,
                                      int                     q,
                                      int                     shift,
                                      int                     np,
                                      int                     p_first )
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
    clargs->p_first = p_first;
    clargs->me      = me;
    clargs->shift   = shift;
    clargs->np_inv  = np;

    rt_starpu_insert_task(
        &cl_zperm_allreduce_col,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_col_args_s),
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, m,  me),
        STARPU_R,                 RTBLKADDR(U, CHAMELEON_Complex64_t, m,  src),
        STARPU_R,                 ipiv_handle,
        STARPU_EXECUTE_ON_NODE,   me,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, m, src) );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
INSERT_TASK_zperm_allreduce_send_col( const RUNTIME_option_t *options,
                                      CHAM_desc_t            *U,
                                      int                     me,
                                      int                     dst,
                                      int                     m )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zperm_allreduce_send_col, 1 );

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, m, me ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
}

static void
INSERT_TASK_zperm_allreduce_recv_col( const RUNTIME_option_t *options,
                                      cham_dir_t              dir,
                                      CHAM_desc_t            *U,
                                      CHAM_ipiv_t            *ipiv,
                                      int                     ipivk,
                                      int                     me,
                                      int                     src,
                                      int                     m,
                                      int                     tempnn,
                                      int                     p,
                                      int                     q,
                                      int                     shift,
                                      int                     np,
                                      int                     p_first )
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

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_allreduce_send_col, zperm_allreduce_col, redux_col, 3 );

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
    clargs->p_first = p_first;
    clargs->me      = me;
    clargs->shift   = shift;
    clargs->np_inv  = np;

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
        chameleon_error( "INSERT_TASK_zperm_allreduce", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, m, src) );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
zperm_allreduce_col_chameleon_starpu_task( const RUNTIME_option_t     *options,
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
    int  np_iter       = np_involved;
    int  p_recv, p_send, me, p_first;
    int  shift = 1;

    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
    }
    else {
        p_first = proc_involved[0];
        for( me = 0; me < np_involved; me++ ) {
            if ( proc_involved[me] == A->myrank ) {
                break;
            }
        }
        assert( me < np_involved );
        while ( np_iter > 1 ) {
            p_send = proc_involved[ ( me + shift               ) % np_involved ];
            p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ];

            INSERT_TASK_zperm_allreduce_send_col( options, U, A->myrank, p_send, m );
            INSERT_TASK_zperm_allreduce_recv_col( options, dir, U, ipiv, ipivk, A->myrank, p_recv,
                                                  m, k == (A->nt-1) ? A->n - k * A->nb : A->nb,
                                                  chameleon_desc_datadist_get_iparam(A, 0),
                                                  chameleon_desc_datadist_get_iparam(A, 1),
                                                  shift, np_involved, p_first );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}

void
INSERT_TASK_zperm_allreduce_col( const RUNTIME_option_t *options,
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
        zperm_allreduce_col_chameleon_starpu_task( options, dir, A, U, Um, Un, ipiv, ipivk, m, k, &(tmp->reduce) );
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

#else

void
INSERT_TASK_zperm_allreduce_send_invp_col( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     m,
                                           int                     k )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)m;
    (void)k;
}

void
INSERT_TASK_zperm_allreduce_col( const RUNTIME_option_t *options,
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
    (void)m;
    (void)k;
    (void)ws;
}

#endif
