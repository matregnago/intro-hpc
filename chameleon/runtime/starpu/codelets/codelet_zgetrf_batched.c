/**
 *
 * @file starpu/codelet_zgetrf_batched.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel batched StarPU codelets
 *
 * @version 1.4.0
 * @comment Codelets to perform batched panel factorization with partial pivoting
 *
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgetrf_batched_args_s {
    const char              *cl_name;
    int                      tasks_nbr;
    int                      diag;
    int                      readUp;
    int                      h;
    int                      ib;
    int                      m[CHAMELEON_BATCH_SIZE];
    int                      n[CHAMELEON_BATCH_SIZE];
    int                      m0[CHAMELEON_BATCH_SIZE];
    struct starpu_data_descr handle_mode[CHAMELEON_BATCH_SIZE];
};

static inline double flops_zgetrf_percol_batched( int *m, int *n, int h, int t )
{
    double flops = 0.;
    int k;
    for ( k = 0; k < t; k ++ ) {
        flops += flops_zgetrf_percol_offdiag( m[k], n[k], h );
    }
    return flops;
}

static inline double flops_zgetrf_blocked_batched( int *m, int *n, int h, int ib, int d, int t )
{
    double flops = 0.;
    int k;
    if ( d == 1 ) {
        flops += flops_zgetrf_blocked_diag( m[0]-h, n[0], h, ib );
    }
    for ( k = d; k < t; k ++ ) {
        flops += flops_zgetrf_blocked_offdiag( m[k]-h, n[k], h, ib );
    }
    return flops;
}

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgetrf_panel_offdiag_batched_cpu_func( void *descr[],
                                          void *cl_arg )
{
    struct cl_zgetrf_batched_args_s *clargs  = (struct cl_zgetrf_batched_args_s *) cl_arg;
    cppi_interface_t                *nextpiv = (cppi_interface_t*) descr[ clargs->tasks_nbr ];
    cppi_interface_t                *prevpiv = (cppi_interface_t*) descr[ clargs->tasks_nbr + 1 ];
    int                              i, m, n, h, m0, lda;
    CHAM_tile_t                     *tileA;

    nextpiv->h = clargs->h;
    nextpiv->has_diag = chameleon_max( -1, nextpiv->has_diag );

    for ( i = 0; i < clargs->tasks_nbr; i++ ) {
        tileA = cti_interface_get( descr[ i ] );
        lda   = tileA->ld;
        m     = clargs->m[ i ];
        n     = clargs->n[ i ];
        h     = clargs->h;
        m0    = clargs->m0[ i ];
        CORE_zgetrf_panel_offdiag( m, n, h, m0, n, CHAM_tile_get_ptr(tileA), lda,
                                   NULL, -1, &( nextpiv->pivot ), &( prevpiv->pivot ) );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

CODELETS_CPU( zgetrf_panel_offdiag_batched, cl_zgetrf_panel_offdiag_batched_cpu_func )

void
INSERT_TASK_zgetrf_panel_offdiag_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          void **clargs_ptr,
                                          CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_offdiag_batched: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    int    task_num   = 0;
    int    batch_size = ((struct chameleon_pzgetrf_s *)ws)->batch_size_blas2;
    struct cl_zgetrf_batched_args_s *clargs = *clargs_ptr;
    int rankA = A->get_rankof( A, Am, An );
    if ( rankA != A->myrank ) {
        return;
    }

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zgetrf_batched_args_s ) ) ;
        memset( clargs, 0, sizeof( struct cl_zgetrf_batched_args_s ) );
        clargs->tasks_nbr   = 0;
        clargs->h           = h;
        clargs->cl_name     = "zgetrf_panel_offdiag_batched";

        *clargs_ptr = clargs;
    }

    task_num               = clargs->tasks_nbr;
    clargs->m[ task_num ]  = m;
    clargs->n[ task_num ]  = n;
    clargs->m0[ task_num ] = m0;
    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;
    /* Refine name */
    clargs->cl_name = chameleon_codelet_name( clargs->cl_name, 1,
                                              A->get_blktile( A, Am, An ) );

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zgetrf_panel_offdiag_batched_flush( options, A, An, clargs_ptr, pivot );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void
INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                void **clargs_ptr,
                                                CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_offdiag_batched: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    struct cl_zgetrf_batched_args_s *clargs = *clargs_ptr;
    int rankA = A->myrank;

    if ( clargs == NULL ) {
        return;
    }
    int access_npiv = ( clargs->h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( clargs->h == 0 )        ? STARPU_NONE : STARPU_R;

    rt_starpu_insert_task(
        &cl_zgetrf_panel_offdiag_batched,
        /* Task codelet arguments */
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zgetrf_batched_args_s),
        STARPU_DATA_MODE_ARRAY,   clargs->handle_mode, clargs->tasks_nbr,
        access_npiv,              RUNTIME_pivot_getaddr( pivot, rankA, clargs->h   ),
        access_ppiv,              RUNTIME_pivot_getaddr( pivot, rankA, clargs->h-1 ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        STARPU_NAME,              clargs->cl_name,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)An;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void
INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                void **clargs_ptr,
                                                CHAM_desc_pivot_t *pivot )
{
    struct cl_zgetrf_batched_args_s *myclargs = *clargs_ptr;
    int rankA = A->myrank;
    int k, ret, access_npiv, access_ppiv;
    struct starpu_task *task;

    if ( myclargs == NULL ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_panel_percol_offdiag_batched, zgetrf_panel_offdiag_batched, zgetrf_batched, myclargs->tasks_nbr + 2 );

    access_npiv = ( myclargs->h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    access_ppiv = ( myclargs->h == 0 )        ? STARPU_NONE : STARPU_R;

    /*
     * Register the data handles, no exchange needed
     */
    starpu_cham_exchange_init_params( options, &params, rankA );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_RW );
    }
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, myclargs->h ),   access_npiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, myclargs->h-1 ), access_ppiv );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = flops_zgetrf_percol_batched( myclargs->m, myclargs->n, myclargs->h, myclargs->tasks_nbr );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_percol_diag", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
    (void)An;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgetrf_panel_blocked_batched_cpu_func( void *descr[],
                                          void *cl_arg )
{
    struct cl_zgetrf_batched_args_s *clargs = ( struct cl_zgetrf_batched_args_s * ) cl_arg;
    int                            *ipiv;
    cppi_interface_t               *nextpiv = (cppi_interface_t*) descr[clargs->tasks_nbr ];
    cppi_interface_t               *prevpiv = (cppi_interface_t*) descr[clargs->tasks_nbr + 1];
    int                             i, h, ib;
    CHAM_tile_t                    *tileA, *tileU;
    CHAMELEON_Complex64_t          *U   = NULL;
    int                             ldu = -1;

    nextpiv->h = clargs->h;
    nextpiv->has_diag = chameleon_max( -1, nextpiv->has_diag);

    h  = clargs->h;
    ib = clargs->ib;
    i  = 0;
    if ( clargs->diag ) {
        if ( h == 0 ) {
            ipiv = (int *)STARPU_VECTOR_GET_PTR(descr[clargs->tasks_nbr + 1]);
        }
        else {
            ipiv = (int *)STARPU_VECTOR_GET_PTR(descr[clargs->tasks_nbr + 2]);
        }
        if ( clargs->readUp ) {
            tileU = cti_interface_get( descr[ clargs->tasks_nbr + 3 ] );
            U     = CHAM_tile_get_ptr( tileU );
            ldu   = tileU->ld;
        }
        tileA             = cti_interface_get( descr[ 0 ] );
        nextpiv->has_diag = 1;
        CORE_zgetrf_panel_diag( clargs->m[i], clargs->n[i], h, clargs->m0[i], ib,
                                CHAM_tile_get_ptr( tileA ), tileA->ld,
                                U, ldu,
                                ipiv, &(nextpiv->pivot), &(prevpiv->pivot) );
        i++;
    }
    if ( clargs->readUp ) {
        tileU = cti_interface_get( descr[ clargs->tasks_nbr + 2 + clargs->diag ] );
        U     = CHAM_tile_get_ptr( tileU );
        ldu   = tileU->ld;
    }
    else {
        U   = NULL;
        ldu = -1;
    }
    for ( ; i < clargs->tasks_nbr; i++ ) {
        tileA = cti_interface_get( descr[ i ] );
        CORE_zgetrf_panel_offdiag( clargs->m[i], clargs->n[i], h, clargs->m0[i], ib,
                                   CHAM_tile_get_ptr(tileA), tileA->ld,
                                   U, ldu,
                                   &( nextpiv->pivot ), &( prevpiv->pivot ) );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

CODELETS_CPU( zgetrf_panel_blocked_batched, cl_zgetrf_panel_blocked_batched_cpu_func )

void
INSERT_TASK_zgetrf_panel_blocked_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0, int readUp,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          CHAM_desc_t *U, int Um, int Un,
                                          void **clargs_ptr,
                                          CHAM_ipiv_t *ipiv,
                                          CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_panel_blocked_batched: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    struct chameleon_pzgetrf_s *tmp = (struct chameleon_pzgetrf_s *) ws;
    int ib         = tmp->ib;
    int batch_size = tmp->batch_size;
    int task_num   = 0;
    struct cl_zgetrf_batched_args_s *clargs = *clargs_ptr;

#if defined ( CHAMELEON_USE_MPI )
    int rankA = A->get_rankof(A, Am, An);
    if ( rankA != A->myrank ) {
        return;
    }
#endif

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zgetrf_batched_args_s ) );
        memset( clargs, 0, sizeof( struct cl_zgetrf_batched_args_s ) );
        clargs->tasks_nbr         = 0;
        clargs->diag              = ( Am == An );
        clargs->readUp            = readUp;
        clargs->ib                = ib;
        clargs->h                 = h;
        clargs->cl_name           = "zgetrf_panel_blocked_batched";

        *clargs_ptr = clargs;
    }

    task_num                               = clargs->tasks_nbr;
    clargs->m[ task_num ]                  = m;
    clargs->n[ task_num ]                  = n;
    clargs->m0[ task_num ]                 = m0;
    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;
    /* Refine name */
    clargs->cl_name = chameleon_codelet_name( clargs->cl_name, 1,
                                              A->get_blktile( A, Am, An ) );

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zgetrf_panel_blocked_batched_flush( options, A, An, readUp, U, Um, Un, clargs_ptr, ipiv, pivot );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void
INSERT_TASK_zgetrf_panel_blocked_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An, int readUp,
                                                CHAM_desc_t *U, int Um, int Un,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv,
                                                CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_panel_blocked_batched: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    int accessU, access_npiv, access_ipiv, access_ppiv;
    struct cl_zgetrf_batched_args_s *clargs = *clargs_ptr;
    int rankA = A->myrank;
    assert( rankA == Um );

    if ( clargs == NULL ) {
        return;
    }

    access_npiv = ( clargs->h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    access_ipiv = ( clargs->h == 0 )        ? STARPU_W    : STARPU_RW;
    access_ppiv = ( clargs->h == 0 )        ? STARPU_NONE : STARPU_R;
    accessU     = ( readUp )                ? STARPU_R    : STARPU_NONE;

    /* If there isn't a diag task then use offdiag access */
    if ( clargs->diag == 0 ) {
        access_ipiv = STARPU_NONE;
    }

    rt_starpu_insert_task(
        &cl_zgetrf_panel_blocked_batched,
        /* Task codelet arguments */
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zgetrf_batched_args_s),
        STARPU_DATA_MODE_ARRAY,   clargs->handle_mode, clargs->tasks_nbr,
        access_npiv,              RUNTIME_pivot_getaddr( pivot, rankA, clargs->h ),
        access_ppiv,              RUNTIME_pivot_getaddr( pivot, rankA, clargs->h - 1 ),
        access_ipiv,              RUNTIME_ipiv_getaddr( ipiv, An ),
        accessU,                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        STARPU_NAME,              clargs->cl_name,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void
INSERT_TASK_zgetrf_panel_blocked_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An, int readUp,
                                                CHAM_desc_t *U, int Um, int Un,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv,
                                                CHAM_desc_pivot_t *pivot )
{
    struct cl_zgetrf_batched_args_s *myclargs = *clargs_ptr;
    int rankA = A->myrank;
    int accessU, access_npiv, access_ipiv, access_ppiv, k;
    int ret;
    struct starpu_task *task;

    if ( myclargs == NULL ) {
        return;
    }
    assert( rankA == Um );

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_panel_blocked_batched, zgetrf_panel_blocked_batched, zgetrf_batched, myclargs->tasks_nbr + 4 );

    access_npiv = ( myclargs->h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    access_ipiv = ( myclargs->h == 0 )        ? STARPU_W    : STARPU_RW;
    access_ppiv = ( myclargs->h == 0 )        ? STARPU_NONE : STARPU_R;
    accessU     = ( readUp )                  ? STARPU_R    : STARPU_NONE;

    /* If there isn't a diag task then use offdiag access */
    if ( myclargs->diag == 0 ) {
        access_ipiv = STARPU_NONE;
    }

    /*
     * Register the data handles, exchange needed only for U
     */
    starpu_cham_exchange_init_params( options, &params, rankA );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_RW );
    }
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, myclargs->h ),   access_npiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, myclargs->h-1 ), access_ppiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_ipiv_getaddr( ipiv, An),                          access_ipiv );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),              accessU );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = flops_zgetrf_blocked_batched( myclargs->m, myclargs->n, myclargs->h, myclargs->ib,
                                                myclargs->diag, myclargs->tasks_nbr );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_panel_blocked_batched", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
