/**
 *
 * @file starpu/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
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

struct cl_zlaswp_single_tile_s {
    int m0;
    int m;
    int n;
    int perm_m;
};

struct cl_zlaswp_batched_args_s {
    int                            side;
    int                            tasks_nbr;
    int                            k;
    int                            perm_mt;
    struct cl_zlaswp_single_tile_s tiles[CHAMELEON_BATCH_SIZE];
    struct starpu_data_descr       handle_mode[CHAMELEON_BATCH_SIZE];
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaswp_batched_cpu_func( void *descr[],
                            void *cl_arg )
{
    int          i, *permget, *permset;
    CHAM_tile_t *A, *WAP, *WA;
    struct cl_zlaswp_batched_args_s *clargs  = ( struct cl_zlaswp_batched_args_s * ) cl_arg;
    struct cl_zlaswp_single_tile_s  *tilearg = clargs->tiles;
    starpu_cham_tile_interface_t   **descrA;

    permget = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    permset = (int *)STARPU_VECTOR_GET_PTR( descr[1] );
    WAP     = (CHAM_tile_t *)cti_interface_get( descr[2] );
    WA      = (CHAM_tile_t *)cti_interface_get( descr[3] );

    descrA = (starpu_cham_tile_interface_t**)&(descr[4]);
    if ( clargs->perm_mt < 0 ) {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_get( clargs->side, tilearg->m0, tilearg->m,
                              tilearg->n, clargs->k, A, WAP, permget );
            TCORE_zlaswp_set( clargs->side, tilearg->m0, tilearg->m,
                              tilearg->n, clargs->k, WA, A, permset );

            descrA++;
            tilearg++;
        }
    }
    else {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_get_idx( clargs->side, tilearg->m0, tilearg->m,
                                  tilearg->n, clargs->k, A, WAP,
                                  tilearg->perm_m, clargs->perm_mt, permget );
            TCORE_zlaswp_set_idx( clargs->side, tilearg->m0, tilearg->m,
                                  tilearg->n, clargs->k, WA, A,
                                  tilearg->perm_m, clargs->perm_mt, permset );

            descrA++;
            tilearg++;
        }
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_batched, cl_zlaswp_batched_cpu_func )

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 cham_side_t             side,
                                 cham_dir_t              dir,
                                 int                     m0,
                                 int                     m,
                                 int                     n,
                                 int                     k,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv, int ipivk,
                                 const CHAM_desc_t      *A,   int Am,   int An,
                                 const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                 const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                 void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzlaswp_s *)ws)->batch_size_swap;
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    if( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zlaswp_batched_args_s ) ) ;
        clargs->side      = side;
        clargs->tasks_nbr = 0;
        clargs->k         = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_mt = A->mt - ipivk;
            }
            else {
                clargs->perm_mt = A->nt - ipivk;
            }
        }
        else {
            clargs->perm_mt = -1;
        }
        *clargs_ptr = clargs;
    }

    task_num = clargs->tasks_nbr;
    clargs->tiles[ task_num ].m0 = m0;
    clargs->tiles[ task_num ].m  = m;
    clargs->tiles[ task_num ].n  = n;
    if ( ipiv->withidx ) {
        if ( side == ChamLeft ) {
            clargs->tiles[ task_num ].perm_m = Am - ipivk;
        }
        else {
            clargs->tiles[ task_num ].perm_m = An - ipivk;
        }
    }
    else {
        clargs->tiles[ task_num ].perm_m = -1;
    }

    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zlaswp_batched_flush( options, dir, ipiv, ipivk, WA, WAm, WAn, WAP, WAPm, WAPn, clargs_ptr );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                       const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                       void                  **clargs_ptr )
{
    struct cl_zlaswp_batched_args_s *clargs   = *clargs_ptr;
    int                              nhandles;
    void                            *ipiv_handle_get;
    void                            *ipiv_handle_set;

    if( clargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle_get = RUNTIME_ipiv_getperm( ipiv, ipivk );
        ipiv_handle_set = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle_get = RUNTIME_ipiv_getinvp( ipiv, ipivk );
        ipiv_handle_set = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    nhandles = clargs->tasks_nbr;
    rt_starpu_insert_task(
        &cl_zlaswp_batched,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswp_batched_args_s),
        STARPU_R,                   ipiv_handle_get,
        STARPU_R,                   ipiv_handle_set,
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(WAP, ChamComplexDouble, WAPm, WAPn),
        STARPU_R,                   RTBLKADDR(WA,  ChamComplexDouble, WAm,  WAn),
        STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *WA,   int WAm,  int WAn,
                                       const CHAM_desc_t      *WAP,  int WAPm, int WAPn,
                                       void                  **clargs_ptr )
{
    int                              ret, k;
    struct starpu_task              *task;
    struct cl_zlaswp_batched_args_s *myclargs = *clargs_ptr;
    void                            *ipiv_handle_get;
    void                            *ipiv_handle_set;

    if( myclargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle_get = RUNTIME_ipiv_getperm( ipiv, ipivk );
        ipiv_handle_set = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle_get = RUNTIME_ipiv_getinvp( ipiv, ipivk );
        ipiv_handle_set = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_batched, myclargs->tasks_nbr + 4 );

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, WA->myrank );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle_get,
                                                  STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle_set,
                                                  STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( WAP, ChamComplexDouble, WAPm, WAPn ),
                                                  STARPU_RW | STARPU_COMMUTE );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( WA, ChamComplexDouble, WAm, WAn ),
                                                  STARPU_R );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_RW );
    }

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_batched", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaswp_get_batched_cpu_func( void *descr[],
                                void *cl_arg )
{
    int          i, *perm;
    CHAM_tile_t *A, *WAP;
    struct cl_zlaswp_batched_args_s *clargs  = ( struct cl_zlaswp_batched_args_s * ) cl_arg;
    struct cl_zlaswp_single_tile_s  *tilearg = clargs->tiles;
    starpu_cham_tile_interface_t   **descrA;

    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    WAP  = (CHAM_tile_t *)cti_interface_get( descr[1] );

    descrA = (starpu_cham_tile_interface_t**)&(descr[2]);
    if ( clargs->perm_mt < 0 ) {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_get( clargs->side, tilearg->m0, tilearg->m,
                              tilearg->n, clargs->k, A, WAP, perm );

            descrA++;
            tilearg++;
        }
    }
    else {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_get_idx( clargs->side, tilearg->m0, tilearg->m,
                                  tilearg->n, clargs->k, A, WAP,
                                  tilearg->perm_m, clargs->perm_mt, perm );

            descrA++;
            tilearg++;
        }
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_get_batched, cl_zlaswp_get_batched_cpu_func )

void INSERT_TASK_zlaswp_get_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *A,   int Am,   int An,
                                     const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                     void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzlaswp_s *)ws)->batch_size_swap_get;
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    if ( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zlaswp_batched_args_s ) );
        clargs->side      = side;
        clargs->tasks_nbr = 0;
        clargs->k         = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_mt = A->mt - ipivk;
            }
            else {
                clargs->perm_mt = A->nt - ipivk;
            }
        }
        else {
            clargs->perm_mt = -1;
        }
        *clargs_ptr = clargs;
    }

    task_num = clargs->tasks_nbr;
    clargs->tiles[ task_num ].m0 = m0;
    clargs->tiles[ task_num ].m  = m;
    clargs->tiles[ task_num ].n  = n;
    if ( ipiv->withidx ) {
        if ( side == ChamLeft ) {
            clargs->tiles[ task_num ].perm_m = Am - ipivk;
        }
        else {
            clargs->tiles[ task_num ].perm_m = An - ipivk;
        }
    }
    else {
        clargs->tiles[ task_num ].perm_m = -1;
    }

    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_R;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zlaswp_get_batched_flush( options, dir, ipiv, ipivk, WAP, WAPm, WAPn, clargs_ptr );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_get_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WAP,  int WAPm, int WAPn,
                                           void                  **clargs_ptr )
{
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    int                              nhandles;
    void                            *ipiv_handle;

    if( clargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    nhandles = clargs->tasks_nbr;
    rt_starpu_insert_task(
        &cl_zlaswp_get_batched,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswp_batched_args_s),
        STARPU_R,                   ipiv_handle,
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(WAP, ChamComplexDouble, WAPm, WAPn),
        STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_get_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WAP,  int WAPm, int WAPn,
                                           void                  **clargs_ptr )
{
    int                              ret, k;
    struct starpu_task              *task;
    struct cl_zlaswp_batched_args_s *myclargs = *clargs_ptr;
    void                            *ipiv_handle;

    if( myclargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_get_batched, myclargs->tasks_nbr + 2 );

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, WAP->myrank );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle,
                                                  STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( WAP, ChamComplexDouble, WAPm, WAPn ),
                                                  STARPU_RW | STARPU_COMMUTE );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_R );
    }

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_get_batched", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaswp_set_batched_cpu_func( void *descr[],
                                void *cl_arg )
{
    int          i, *invp;
    CHAM_tile_t *WA, *A;
    struct cl_zlaswp_batched_args_s *clargs  = ( struct cl_zlaswp_batched_args_s * ) cl_arg;
    struct cl_zlaswp_single_tile_s  *tilearg = clargs->tiles;
    starpu_cham_tile_interface_t   **descrA;

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    WA   = (CHAM_tile_t *)cti_interface_get( descr[1] );

    descrA = (starpu_cham_tile_interface_t**)&(descr[2]);
    if ( clargs->perm_mt < 0 ) {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_set( clargs->side, tilearg->m0, tilearg->m,
                              tilearg->n, clargs->k, WA, A, invp );

            descrA++;
            tilearg++;
        }
    }
    else {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrA );

            TCORE_zlaswp_set_idx( clargs->side, tilearg->m0, tilearg->m,
                                  tilearg->n, clargs->k, WA, A,
                                  tilearg->perm_m, clargs->perm_mt, invp );

            descrA++;
            tilearg++;
        }
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_set_batched, cl_zlaswp_set_batched_cpu_func )

void INSERT_TASK_zlaswp_set_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *WA,   int WAm, int WAn,
                                     const CHAM_desc_t      *A,    int Am,  int An,
                                     void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzlaswp_s *)ws)->batch_size_swap_set;
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    if( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zlaswp_batched_args_s ) );
        clargs->side      = side;
        clargs->tasks_nbr = 0;
        clargs->k         = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_mt = A->mt - ipivk;
            }
            else {
                clargs->perm_mt = A->nt - ipivk;
            }
        }
        else {
            clargs->perm_mt = -1;
        }
        *clargs_ptr = clargs;
    }

    task_num = clargs->tasks_nbr;
    clargs->tiles[ task_num ].m0 = m0;
    clargs->tiles[ task_num ].m  = m;
    clargs->tiles[ task_num ].n  = n;
    if ( ipiv->withidx ) {
        if ( side == ChamLeft ) {
            clargs->tiles[ task_num ].perm_m = Am - ipivk;
        }
        else {
            clargs->tiles[ task_num ].perm_m = An - ipivk;
        }
    }
    else {
        clargs->tiles[ task_num ].perm_m = -1;
    }

    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zlaswp_set_batched_flush( options, dir, ipiv, ipivk, WA, WAm, WAn, clargs_ptr );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_set_batched_flush( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *WA,   int WAm, int WAn,
                                       void                  **clargs_ptr )
{
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    int                              nhandles;
    void                            *ipiv_handle;

    if( clargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    nhandles = clargs->tasks_nbr;
    rt_starpu_insert_task(
        &cl_zlaswp_set_batched,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswp_batched_args_s),
        STARPU_R,                   ipiv_handle,
        STARPU_R,                   RTBLKADDR(WA, ChamComplexDouble, WAm, WAn),
        STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_set_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WA,   int WAm, int WAn,
                                           void                  **clargs_ptr )
{
    int                              ret, k;
    struct starpu_task              *task;
    struct cl_zlaswp_batched_args_s *myclargs = *clargs_ptr;
    void                            *ipiv_handle;

    if( myclargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_set_batched, myclargs->tasks_nbr + 2 );

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, WA->myrank );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle,
                                                  STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( WA, ChamComplexDouble, WAm, WAn ),
                                                  STARPU_R );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_RW );
    }

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_set_batched", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
