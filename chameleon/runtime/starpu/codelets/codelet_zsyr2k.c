/**
 *
 * @file starpu/codelet_zsyr2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Loris Lucido
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zsyr2k_args_s {
    cham_uplo_t           uplo;
    cham_trans_t          trans;
    int                   n;
    int                   k;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};


#if !defined(CHAMELEON_SIMULATION)
static void cl_zsyr2k_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zsyr2k_args_s *clargs = (struct cl_zsyr2k_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    TCORE_zsyr2k( clargs->uplo, clargs->trans,
                  clargs->n, clargs->k, clargs->alpha,
                  tileA, tileB, clargs->beta, tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zsyr2k_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    struct cl_zsyr2k_args_s *clargs = (struct cl_zsyr2k_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    TCUDA_zsyr2k( clargs->uplo, clargs->trans,
                  clargs->n, clargs->k, (cuDoubleComplex *)&(clargs->alpha),
                  tileA, tileB, (cuDoubleComplex *)&(clargs->beta), tileC, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void cl_zsyr2k_hip_func(void *descr[], void *cl_arg)
{
    hipblasHandle_t handle = starpu_hipblas_get_local_handle();
    struct cl_zsyr2k_args_s *clargs = (struct cl_zsyr2k_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    HIP_zsyr2k( clargs->uplo, clargs->trans,
                clargs->n, clargs->k,
                (hipDoubleComplex*)&(clargs->alpha), tileA->mat, tileA->ld,
                                                         tileB->mat, tileB->ld,
                (hipDoubleComplex*)&(clargs->beta),  tileC->mat, tileC->ld,
                handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zsyr2k, cl_zsyr2k_cpu_func, cl_zsyr2k_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zsyr2k, cl_zsyr2k_cpu_func, cl_zsyr2k_cuda_func, STARPU_CUDA_ASYNC )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zsyr2k_args_s     *clargs  = NULL;
    int                          exec    = 0;
    const char                  *cl_name = "zsyr2k";
    enum starpu_data_access_mode accessC = STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zsyr2k_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Add commute if possible (rectask is not compatible yet) */
    if ( beta == (CHAMELEON_Complex64_t)1. ) {
        accessC |= STARPU_COMMUTE;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    rt_starpu_insert_task(
        &cl_zsyr2k,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zsyr2k_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,      RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        accessC,       RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zsyr2k ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zsyr2k( k, n ),
        0 );

    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zsyr2k, 3 );
    enum starpu_data_access_mode accessC = STARPU_RW;

    /* Add commute if possible (rectask is not compatible yet) */
    if ( beta == (CHAMELEON_Complex64_t)1. ) {
        accessC |= STARPU_COMMUTE;
    }

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, C->get_rankof( C, Cm, Cn ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, B, Bm, Bn, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, C, Cm, Cn, accessC  );

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
        clargs = malloc( sizeof( struct cl_zsyr2k_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zsyr2k_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zsyr2k_callback );

        /* Flops */
        task->flops = flops_zsyr2k( k, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 3,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ),
                                             C->get_blktile( C, Cm, Cn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zsyr2k", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
