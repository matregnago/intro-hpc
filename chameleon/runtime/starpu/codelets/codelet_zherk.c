/**
 *
 * @file starpu/codelet_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Gwenole Lucas
 * @author Loris Lucido
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zherk_args_s {
    cham_uplo_t  uplo;
    cham_trans_t trans;
    int          n;
    int          k;
    double       alpha;
    double       beta;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zherk_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zherk_args_s *clargs  = (struct cl_zherk_args_s *)(t->cl_arg);
    rectask_args_t         *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t       request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzherk( clargs->uplo, clargs->trans,
                      clargs->alpha, rtargs->tiles[0]->mat,
                      clargs->beta,  rtargs->tiles[1]->mat,
                      rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zherk_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    assert( tileA->flttype == ChamComplexDouble );
    assert( tileC->flttype == ChamComplexDouble );

    TCORE_zherk( clargs->uplo, clargs->trans, clargs->n, clargs->k,
                 clargs->alpha, tileA, clargs->beta, tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zherk_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    TCUDA_zherk( clargs->uplo, clargs->trans, clargs->n, clargs->k,
                 &(clargs->alpha), tileA, &(clargs->beta), tileC, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_zherk_hip_func(void *descr[], void *cl_arg)
{
    hipblasHandle_t handle = starpu_hipblas_get_local_handle();
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    HIP_zherk(
        clargs->uplo, clargs->trans, clargs->n, clargs->k,
        &(clargs->alpha),
        tileA->mat, tileA->ld,
        &(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zherk, cl_zherk_cpu_func, cl_zherk_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zherk, cl_zherk_cpu_func, cl_zherk_cuda_func, STARPU_CUDA_ASYNC )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zherk_args_s      *clargs  = NULL;
    int                          exec    = 0;
    const char                  *cl_name = "zherk";
    enum starpu_data_access_mode accessC = STARPU_RW;
    CHAM_tile_t                 *tileA;
    CHAM_tile_t                 *tileC;
    int                          is_rectask = 0;
    rectask_args_t              *rtargs     = NULL;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    tileA = A->get_blktile( A, Am, An );
    tileC = C->get_blktile( C, Cm, Cn );

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( ( tileA->format & CHAMELEON_TILE_DESC ) &&
                   ( tileC->format & CHAMELEON_TILE_DESC ) );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + sizeof(CHAM_tile_t*) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority ;
        rtargs->tiles[0] = tileA;
        rtargs->tiles[1] = tileC;
        cl_name = "zherk_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zherk_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Add commute if possible (rectask is not compatible yet) */
    if ( (beta == 1.) && !is_rectask ) {
        accessC |= STARPU_COMMUTE;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2, tileA, tileC );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zherk,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zherk_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessC,       RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zherk ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zherk( k, n ),

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zherk )
        0 );

    (void)tileA;
    (void)tileC;
    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zherk, 2 );
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
        clargs = malloc( sizeof( struct cl_zherk_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zherk_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zherk_callback );

        /* Flops */
        task->flops = flops_zherk( k, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             C->get_blktile( C, Cm, Cn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zherk", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
