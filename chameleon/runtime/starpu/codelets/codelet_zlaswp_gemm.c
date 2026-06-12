/**
 *
 * @file starpu/codelet_zlaswp_gemm.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp_set and gemm on a panel
 *
 * @version 1.4.0
 * @author Alycia Lisito
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

struct cl_zlaswp_gemm_args_s {
    int                            side;
    int                            tasks_nbr;
    int                            k;
    int                            perm_mt;
    struct cl_zlaswp_single_tile_s tiles[CHAMELEON_BATCH_SIZE];
    struct starpu_data_descr       handle_mode[2*CHAMELEON_BATCH_SIZE];
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaswp_gemm_cpu_func( void *descr[],
                         void *cl_arg )
{
    int          i, *invp;
    CHAM_tile_t *WA, *A, *B, *C;
    struct cl_zlaswp_gemm_args_s   *clargs  = ( struct cl_zlaswp_gemm_args_s * ) cl_arg;
    struct cl_zlaswp_single_tile_s *tilearg = clargs->tiles;
    starpu_cham_tile_interface_t  **descrAC;
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    WA   = (CHAM_tile_t *)cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *)cti_interface_get( descr[2] );

    descrAC = (starpu_cham_tile_interface_t**)&(descr[3]);
    if ( clargs->perm_mt < 0 ) {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrAC );
            descrAC++;
            C = (CHAM_tile_t *) cti_interface_get( *descrAC );
            descrAC++;

            TCORE_zlaswp_set( clargs->side, tilearg->m0, tilearg->m,
                              tilearg->n, clargs->k, WA, C, invp );
            TCORE_zgemm( ChamNoTrans, ChamNoTrans,
                         tilearg->m, tilearg->n, clargs->k,
                         mzone, A, B,
                         zone,  C );
            tilearg++;
        }
    }
    else {
        for ( i = 0; i < clargs->tasks_nbr; i++ ) {
            A = (CHAM_tile_t *) cti_interface_get( *descrAC );
            descrAC++;
            C = (CHAM_tile_t *) cti_interface_get( *descrAC );
            descrAC++;

            TCORE_zlaswp_set_idx( clargs->side, tilearg->m0, tilearg->m,
                                  tilearg->n, clargs->k, WA, C,
                                  tilearg->perm_m, clargs->perm_mt, invp );
            TCORE_zgemm( ChamNoTrans, ChamNoTrans,
                         tilearg->m, tilearg->n, clargs->k,
                         mzone, A, B,
                         zone,  C );

            tilearg++;
        }
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_gemm, cl_zlaswp_gemm_cpu_func )

void INSERT_TASK_zlaswp_gemm( const RUNTIME_option_t *options,
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
                              const CHAM_desc_t      *B,    int Bm,  int Bn,
                              const CHAM_desc_t      *C,    int Cm,  int Cn,
                              void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzlaswp_s *)ws)->batch_size_swap_set;
    struct cl_zlaswp_gemm_args_s *clargs = *clargs_ptr;
    if ( C->get_rankof( C, Cm, Cn) != C->myrank ) {
        return;
    }

    if( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zlaswp_gemm_args_s ) );
        clargs->side      = side;
        clargs->tasks_nbr = 0;
        clargs->k         = k;
        if ( ipiv->withidx ) {
            if ( side == ChamLeft ) {
                clargs->perm_mt = C->mt - ipivk;
            }
            else {
                clargs->perm_mt = C->nt - ipivk;
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
            clargs->tiles[ task_num ].perm_m = Cm - ipivk;
        }
        else {
            clargs->tiles[ task_num ].perm_m = Cn - ipivk;
        }
    }
    else {
        clargs->tiles[ task_num ].perm_m = -1;
    }

    clargs->handle_mode[ 2 * task_num     ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    clargs->handle_mode[ 2 * task_num     ].mode   = STARPU_R;
    clargs->handle_mode[ 2 * task_num + 1 ].handle = RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn);
    clargs->handle_mode[ 2 * task_num + 1 ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zlaswp_gemm_flush( options, dir, ipiv, ipivk, WA, WAm, WAn,
                                                 B, Bm, Bn, clargs_ptr );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_gemm_flush( const RUNTIME_option_t *options,
                                    cham_dir_t              dir,
                                    const CHAM_ipiv_t      *ipiv, int ipivk,
                                    const CHAM_desc_t      *WA,   int WAm, int WAn,
                                    const CHAM_desc_t      *B,    int Bm,  int Bn,
                                    void                  **clargs_ptr )
{
    struct cl_zlaswp_gemm_args_s *clargs = *clargs_ptr;
    int   nhandles;
    void *ipiv_handle;

    if( clargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    nhandles = clargs->tasks_nbr * 2;
    rt_starpu_insert_task(
        &cl_zlaswp_gemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlaswp_gemm_args_s),

        /* Task handles */
        STARPU_R,               ipiv_handle,
        STARPU_R,               RTBLKADDR(WA, ChamComplexDouble, WAm, WAn),
        STARPU_R,               RTBLKADDR(B,  ChamComplexDouble, Bm,  Bn),
        STARPU_DATA_MODE_ARRAY, clargs->handle_mode, nhandles,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_gemm_flush( const RUNTIME_option_t *options,
                                    cham_dir_t              dir,
                                    const CHAM_ipiv_t      *ipiv, int ipivk,
                                    const CHAM_desc_t      *WA,   int WAm, int WAn,
                                    const CHAM_desc_t      *B,    int Bm,  int Bn,
                                    void                  **clargs_ptr )
{
    int                 ret, k;
    struct starpu_task *task;
    struct cl_zlaswp_gemm_args_s *myclargs = *clargs_ptr;
    void *ipiv_handle;

    if( myclargs == NULL ) {
        return;
    }

    if ( dir == ChamDirForward ){
        ipiv_handle = RUNTIME_ipiv_getinvp( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_ipiv_getperm( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_gemm, myclargs->tasks_nbr * 2 + 3 );

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
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( B, ChamComplexDouble, Bm, Bn ),
                                                  STARPU_R );
    for ( k = 0; k < myclargs->tasks_nbr * 2; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, myclargs->handle_mode[ k ].mode );
    }

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_gemm_args_s );
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
