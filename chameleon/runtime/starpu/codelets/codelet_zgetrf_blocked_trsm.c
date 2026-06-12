/**
 *
 * @file starpu/codelet_zgetrf_blocked_trsm.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel StarPU codelets
 *
 * @version 1.4.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgetrf_blocked_args_s {
        int                 m;
        int                 n;
        int                 h;
        int                 m0;
        int                 ib;
        RUNTIME_sequence_t *sequence;
        RUNTIME_request_t  *request;
};

struct cl_zgetrf_cpy_args_s {
    int ib;
    int h;
    int n;
};

CHAMELEON_CL_CB( zgetrf_blocked_trsm, cti_handle_get_m(task->handles[0]), 0, 0, M )

#if !defined(CHAMELEON_SIMULATION)

static void
zgetrf_cpy_pivrow_in_Up( CHAM_tile_t      *tileU,
                         cppi_interface_t *prevpiv,
                         int               h,
                         int               n,
                         int               ib )
{
    CHAMELEON_Complex64_t *U = CHAM_tile_get_ptr( tileU );
    int                    ldu = tileU->ld;

    coreblas_kernel_trace( tileU );

    /* Copy the final max line of the block and solve */
    cblas_zcopy( n, prevpiv->pivot.pivrow, 1,
                    U + ( h % ib ), ldu );
}

static void cl_zgetrf_cpy_pivrow_in_Up_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_cpy_args_s *clargs = (struct cl_zgetrf_cpy_args_s *)cl_arg;
    CHAM_tile_t      *tileU = cti_interface_get(descr[0]);
    cppi_interface_t *prevpiv = (cppi_interface_t*) descr[1];

    zgetrf_cpy_pivrow_in_Up( tileU, prevpiv, clargs->h, clargs->n, clargs-> ib );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

CODELETS_CPU(zgetrf_cpy_pivrow_in_Up, cl_zgetrf_cpy_pivrow_in_Up_cpu_func)

#if defined(CHAMELEON_STARPU_USE_INSERT) /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_cpy_pivrow_in_Up( const RUNTIME_option_t *options,
                                          CHAM_desc_t            *Up,
                                          int                     Upm,
                                          int                     k,
                                          int                     h,
                                          int                     ib,
                                          int                     n,
                                          CHAM_desc_pivot_t      *pivot )
{
    struct cl_zgetrf_cpy_args_s *clargs = malloc( sizeof( struct cl_zgetrf_cpy_args_s ) );
    int accessUp = ( ( ( h % ib ) == 0 ) || ( ib == 1 ) ) ? STARPU_W : STARPU_RW;
    clargs->ib = ib;
    clargs->h  = h;
    clargs->n  = n;

    assert( Up->myrank == Up->get_rankof( Up, Upm, 0 ) );
    assert( Up->myrank == Upm );
    assert( ( ( h + 1 ) % ib ) != 0 );

    rt_starpu_insert_task(
        &cl_zgetrf_cpy_pivrow_in_Up,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zgetrf_cpy_args_s),
        accessUp,                 RTBLKADDR(Up, CHAMELEON_Complex64_t, Upm, 0),
        STARPU_R,                 RUNTIME_pivot_getaddr( pivot, Up->myrank, h ),
        STARPU_EXECUTE_ON_NODE,   Up->myrank,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    (void)k;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_cpy_pivrow_in_Up( const RUNTIME_option_t *options,
                                          CHAM_desc_t            *Up,
                                          int                     Upm,
                                          int                     k,
                                          int                     h,
                                          int                     ib,
                                          int                     n,
                                          CHAM_desc_pivot_t      *pivot )
{
    int    ret;
    struct starpu_task *task;
    int    accessUp = ( ( ( h % ib ) == 0 ) || ( ib == 1 ) ) ? STARPU_W : STARPU_RW;
    int    rankUp   = Up->get_rankof(Up, Upm, 0);
    assert( Up->myrank == Upm );
    assert( Up->myrank == rankUp );
    assert( ( ( h + 1 ) % ib ) != 0 );
    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_cpy_pivrow_in_Up, zgetrf_cpy_pivrow_in_Up, zgetrf_cpy, 2 )

    starpu_cham_exchange_init_params( options, &params, rankUp );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( Up, ChamComplexDouble, Upm, 0 ),     accessUp );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, Up->myrank, h ), STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs     = malloc( sizeof( struct cl_zgetrf_cpy_args_s ) );
    clargs->ib = ib;
    clargs->h  = h;
    clargs->n  = n;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_cpy_args_s );
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
        chameleon_error( "INSERT_TASK_zgetrf_cpy_pivrow_in_Up", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)k;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */


#if !defined(CHAMELEON_SIMULATION)
static const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

static void cl_zgetrf_blocked_trsm_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_blocked_args_s *clargs = (struct cl_zgetrf_blocked_args_s *)cl_arg;
    CHAM_tile_t           *tileU;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U;
    int                    ldu;

    tileU   = cti_interface_get(descr[0]);
    prevpiv = (cppi_interface_t*) descr[1];
    U       = CHAM_tile_get_ptr( tileU );
    ldu     = tileU->ld;

    coreblas_kernel_trace( tileU );

    /* Copy the final max line of the block and solve */
    cblas_zcopy( clargs->n, prevpiv->pivot.pivrow, 1,
                            U +  clargs->m - 1, ldu );

    if ( ( clargs->n - clargs->h ) > 0 ) {
        cblas_ztrsm( CblasColMajor,
                     CblasLeft, CblasLower,
                     CblasNoTrans, CblasUnit,
                     clargs->ib, clargs->n - clargs->h,
                     CBLAS_SADDR(zone), U + (clargs->h-clargs->ib) * ldu, ldu,
                                        U +  clargs->h             * ldu, ldu );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_blocked_trsm, cl_zgetrf_blocked_trsm_cpu_func)

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int                     m,
                                      int                     n,
                                      int                     h,
                                      int                     ib,
                                      CHAM_desc_t            *U,
                                      int                     Um,
                                      int                     Un,
                                      CHAM_desc_pivot_t      *pivot )
{
    const char *cl_name = "zgetrf_blocked_trsm";
    int rankU = U->get_rankof(U, Um, Un);

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( U, Um, Un );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      U->get_blktile( U, Um, Un ) );

    if ( U->myrank != rankU ) {
        return;
    }

    /* Set codelet parameters */
    struct cl_zgetrf_blocked_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m  = m;
    clargs->n  = n;
    clargs->h  = h;
    clargs->ib = ib;

    rt_starpu_insert_task(
        &cl_zgetrf_blocked_trsm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgetrf_blocked_args_s),

        /* Task handles */
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_R,                 RUNTIME_pivot_getaddr( pivot, rankU, h-1 ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgetrf_blocked_trsm ),
        STARPU_NAME,              cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int                     m,
                                      int                     n,
                                      int                     h,
                                      int                     ib,
                                      CHAM_desc_t            *U,
                                      int                     Um,
                                      int                     Un,
                                      CHAM_desc_pivot_t      *pivot )
{
    int ret;
    struct starpu_task *task;
    int rankU = U->get_rankof(U, Um, Un);

    if ( U->myrank != rankU ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_blocked_trsm, zgetrf_blocked_trsm, zgetrf_blocked, 2 );

    /*
     * Register the data handles, no exchange needed
     */
    starpu_cham_exchange_init_params( options, &params, rankU );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( U, ChamComplexDouble, Um, Un ),      STARPU_RW );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankU, h-1 ), STARPU_R  );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m  = m;
    clargs->n  = n;
    clargs->h  = h;
    clargs->ib = ib;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_blocked_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgetrf_blocked_trsm_callback );

    /* Flops */
    task->flops = flops_zgetrf_trsm( m, n, h, ib );

    /* Refine name */
    task->name = chameleon_codelet_name( cl_name, 1, U->get_blktile( U, Um, Un ) );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_blocked_diag", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
