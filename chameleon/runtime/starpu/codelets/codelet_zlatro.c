/**
 *
 * @file starpu/codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlatro_args_s {
    cham_uplo_t  uplo;
    cham_trans_t trans;
    int          m;
    int          n;
};

static cham_fixdbl_t
flops_zlatro( cham_uplo_t uplo, int _M, int _N )
{
    cham_fixdbl_t  flops = 0.;
    cham_fixdbl_t  M     = _M;
    cham_fixdbl_t  N     = _N;

    switch( uplo ) {
        case ChamUpperLower:
            flops = M * N;
            break;
        case ChamLower:
            flops = N * (M+1) * M / 2;
            break;
        case ChamUpper:
            flops = M * (N+1) * N / 2;
            break;
    }
    return flops;
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlatro_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlatro_args_s *clargs = (struct cl_zlatro_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_zlatro( clargs->uplo, clargs->trans, clargs->m, clargs->n, tileA, tileB );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zlatro_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_zlatro_args_s *clargs = (struct cl_zlatro_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    cublasHandle_t handle = starpu_cublas_get_local_handle();

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCUDA_zlatro( clargs->uplo, clargs->trans, clargs->m, clargs->n, tileA, tileB, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zlatro, cl_zlatro_cpu_func, cl_zlatro_cuda_func, STARPU_CUDA_ASYNC )

#if defined(CHAMELEON_STARPU_USE_INSERT)
/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_zlatro_args_s *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "zlatro";


    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlatro_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
    }


    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ) );

    rt_starpu_insert_task(
        &cl_zlatro,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zlatro_args_s),
        STARPU_R,                 RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,                 RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zlatro ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zlatro( uplo, m, n ),
        0 );

    (void)mb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    int                 ret;
    struct starpu_task *task;

    INSERT_TASK_COMMON_PARAMETERS( zlatro, 2 );

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, B, Bm, Bn, STARPU_W );
    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        return;
    }

    if ( params.do_execute )
    {
        task = starpu_task_create();
        task->cl = cl;

        clargs = malloc( sizeof( struct cl_zlatro_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zlatro_args_s );
        task->cl_arg_free = 1;

        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zlatro_callback );

        /* Flops */
        task->flops = flops_zlatro( uplo, m, n );

        /* Refine name */
        task->name = cl_name;

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zlatro", "Failed to submit the task to StarPU" );
            return;
        }
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)mb;
}

#endif
