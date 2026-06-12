/**
 *
 * @file starpu/codelet_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgeadd_args_s {
    cham_trans_t          trans;
    int                   m;
    int                   n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgeadd_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zgeadd_args_s *clargs = (struct cl_zgeadd_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_zgeadd( clargs->trans, clargs->m, clargs->n,
                  clargs->alpha, tileA, clargs->beta, tileB );

    return;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zgeadd_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_zgeadd_args_s *clargs = (struct cl_zgeadd_args_s *)cl_arg;
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    CHAM_tile_t   *tileA;
    CHAM_tile_t   *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );

    TCUDA_zgeadd( clargs->trans, clargs->m, clargs->n,
                  (cuDoubleComplex *)&(clargs->alpha), tileA, (cuDoubleComplex *)&(clargs->beta), tileB, handle );

    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgeadd, cl_zgeadd_cpu_func, cl_zgeadd_cuda_func, STARPU_CUDA_ASYNC )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    struct cl_zgeadd_args_s *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "zgeadd";
    int                      accessB;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zgeadd_args_s ) );
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Reduce the B access if needed */
    accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgeadd,

        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgeadd_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,  RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgeadd ),
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zgeadd, 2 );
    int accessB;

    /* Reduce the B access if needed */
    accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, B, Bm, Bn, accessB  );

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        return;
    }

    if ( params.do_execute )
    {
        int ret;
        struct starpu_task *task = starpu_task_create();
        task->cl = cl;

        /* Set codelet parameters */
        clargs = malloc( sizeof( struct cl_zgeadd_args_s ) );
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zgeadd_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgeadd_callback );

        /* Flops */
        //task->flops = flops_zgeadd( m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zgeadd", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
