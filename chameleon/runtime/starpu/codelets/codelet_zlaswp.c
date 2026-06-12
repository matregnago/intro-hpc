/**
 *
 * @file starpu/codelet_zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaswp_args_s {
    cham_side_t side;
    int         m0;
    int         m;
    int         n;
    int         k;
    int         perm_m;
    int         perm_mt;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_get_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswp_args_s *clargs = (struct cl_zlaswp_args_s *)cl_arg;
    int         *perm;
    CHAM_tile_t *A, *WAP;

    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    WAP  = (CHAM_tile_t *) cti_interface_get( descr[2] );

    if ( clargs->perm_mt < 0 ) {
        TCORE_zlaswp_get( clargs->side, clargs->m0, clargs->m,
                          clargs->n, clargs->k, A, WAP, perm );
    }
    else {
        TCORE_zlaswp_get_idx( clargs->side, clargs->m0, clargs->m,
                              clargs->n, clargs->k, A, WAP,
                              clargs->perm_m, clargs->perm_mt, perm );
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_get, cl_zlaswp_get_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn )
{
    struct cl_zlaswp_args_s *clargs  = NULL;
    const char              *cl_name = "zlaswp_get";
    void                    *ipiv_handle;

    /* Handle cache */
    if ( A->get_rankof( A, Am, An ) != A->myrank ) {
        return;
    }

    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;
    if ( ipiv->withidx ) {
        if ( side == ChamLeft ) {
            clargs->perm_m = Am - ipivk;
        }
        else {
            clargs->perm_m = An - ipivk;
        }
        clargs->perm_mt = ipiv->max_mt - ipivk;
    }
    else {
        clargs->perm_m  = -1;
        clargs->perm_mt = -1;
    }

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      A->get_blktile( A, Am, An ),
                                      WAP->get_blktile( WAP, WAPm, WAPn ) );

    rt_starpu_insert_task(
        &cl_zlaswp_get,

        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlaswp_args_s),

        /* Task handles */
        STARPU_R,                   ipiv_handle,
        STARPU_R,                   RTBLKADDR(A,   ChamComplexDouble, Am,   An  ),
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(WAP, ChamComplexDouble, WAPm, WAPn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        STARPU_NAME, cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn )
{
    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswp_get, zlaswp_get, zlaswp, 3 );
    void *ipiv_handle;

    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, WAP->get_rankof( WAP, WAPm, WAPn ) );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( WAP, ChamComplexDouble, WAPm, WAPn ),
                                STARPU_RW | STARPU_COMMUTE );

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        assert( 0 ); /* We should never end up in tis case */
        return;
    }

    if ( params.do_execute )
    {
        int                 ret;
        struct starpu_task *task = starpu_task_create();
        task->cl = cl;

        /* Set codelet parameters */
        clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
        clargs->side = side;
        clargs->m0   = m0;
        clargs->m    = m;
        clargs->n    = n;
        clargs->k    = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_m = Am - ipivk;
            }
            else {
                clargs->perm_m = An - ipivk;
            }
            clargs->perm_mt = ipiv->max_mt - ipivk;
        }
        else {
            clargs->perm_m  = -1;
            clargs->perm_mt = -1;
        }

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zlaswp_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

        /* Flops */
        task->flops = 0.;

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             WAP->get_blktile( WAP, WAPm, WAPn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zlaswp_get", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_set_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswp_args_s *clargs = (struct cl_zlaswp_args_s *)cl_arg;
    int         *invp;
    CHAM_tile_t *WA, *A;

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    WA   = (CHAM_tile_t *) cti_interface_get( descr[1] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    if ( clargs->perm_mt < 0 ) {
        TCORE_zlaswp_set( clargs->side, clargs->m0, clargs->m,
                          clargs->n, clargs->k, WA, A, invp );
    }
    else {
        TCORE_zlaswp_set_idx( clargs->side, clargs->m0, clargs->m,
                              clargs->n, clargs->k, WA, A,
                              clargs->perm_m, clargs->perm_mt, invp );
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_set, cl_zlaswp_set_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t dir, int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *WA,   int WAm, int WAn,
                             const CHAM_desc_t *A,    int Am,  int An )
{
    struct cl_zlaswp_args_s *clargs  = NULL;
    const char              *cl_name = "zlaswp_set";
    void                    *ipiv_handle;

    /* Handle cache */
    if ( A->get_rankof( A, Am, An) != WA->myrank ) {
        return;
    }

    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;
    if ( ipiv->withidx ) {
        if ( side == ChamLeft ) {
            clargs->perm_m = Am - ipivk;
        }
        else {
            clargs->perm_m = An - ipivk;
        }
        clargs->perm_mt = ipiv->max_mt - ipivk;
    }
    else {
        clargs->perm_m  = -1;
        clargs->perm_mt = -1;
    }

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      WA->get_blktile( WA, WAm, WAn ),
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        &cl_zlaswp_set,
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlaswp_args_s),

        /* Task handles */
        STARPU_R,  ipiv_handle,
        STARPU_R,  RTBLKADDR(WA, ChamComplexDouble, WAm, WAn),
        STARPU_RW, RTBLKADDR(A,  ChamComplexDouble, Am,  An ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        STARPU_NAME, cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *WA, int WAm, int WAn,
                             const CHAM_desc_t *A,  int Am,  int An )
{
    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswp_set, zlaswp_set, zlaswp, 3 );
    void *ipiv_handle;

    if ( A->get_rankof( A, Am, An) != WA->myrank ) {
        return;
    }

    if( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, A->get_rankof( A, Am, An ) );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( WA, ChamComplexDouble, WAm, WAn), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_RW );

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        assert( 0 ); /* We should never end up in tis case */
        return;
    }

    if ( params.do_execute )
    {
        int                 ret;
        struct starpu_task *task = starpu_task_create();
        task->cl = cl;

        /* Set codelet parameters */
        clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
        clargs->side = side;
        clargs->m0   = m0;
        clargs->m    = m;
        clargs->n    = n;
        clargs->k    = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_m = Am - ipivk;
            }
            else {
                clargs->perm_m = An - ipivk;
            }
            clargs->perm_mt = ipiv->max_mt - ipivk;
        }
        else {
            clargs->perm_m  = -1;
            clargs->perm_mt = -1;
        }

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zlaswp_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

        /* Flops */
        task->flops = 0.;

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             WA->get_blktile( WA, WAm, WAn ),
                                             A->get_blktile( A, Am, An ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zlaswp_set", "Failed to submit the task to StarPU" );
            return;
        }
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}
#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if defined(CHAMELEON_USE_MPI)

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_ret_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *A_tile;
    cpui_interface_t      *ws     = (cpui_interface_t*) descr[1];
    cham_side_t            side   = ws->side;
    int                    ldb    = ws->n;
    int                   *index  = ws->ws.index;
    CHAMELEON_Complex64_t *rows   = ws->ws.rows;
    CHAMELEON_Complex64_t *A;
    int                    i, lda, A_inc;

    A_tile = (CHAM_tile_t *) cti_interface_get( descr[0] );
    A      = CHAM_tile_get_ptr( A_tile );

    lda   = ( side == ChamLeft ) ? A_tile->ld : 1;
    A_inc = ( side == ChamLeft ) ? 1          : A_tile->ld;

    for ( i = 0; i < ws->m; i++ ) {
        if ( index[i] != -1 ) {
            cblas_zcopy( ws->n, rows + index[i] * ldb,   1,
                                A    + i        * A_inc, lda );
        }
    }

    (void)cl_arg;
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_ret, cl_zlaswp_ret_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An)
{
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    rt_starpu_insert_task(
        &cl_zlaswp_ret,

        /* Task handles */
        STARPU_W, RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R, RUNTIME_perm_getaddr( ws, Wm, Wn ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An)
{
    int                 rank = A->get_rankof( A, Am, An );
    int                 ret;
    struct starpu_task *task;

    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_ret, 2);

    starpu_cham_exchange_init_params( options, &params, rank );

    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_W );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_perm_getaddr( ws, Wm, Wn ),
                                STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    /* Refine name */
    task->name = cl_name;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_ret", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)clargs;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#endif /* defined(CHAMELEON_USE_MPI) */
