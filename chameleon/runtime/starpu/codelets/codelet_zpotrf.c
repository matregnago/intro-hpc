/**
 *
 * @file starpu/codelet_zpotrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Terry Cojean
 * @author Gwenole Lucas
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zpotrf_args_s {
    cham_uplo_t         uplo;
    int                 n;
    int                 iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zpotrf_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zpotrf_args_s *clargs  = (struct cl_zpotrf_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzpotrf( clargs->uplo, rtargs->tiles[0]->mat,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zpotrf_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zpotrf_args_s *clargs = (struct cl_zpotrf_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    assert( tileA->flttype == ChamComplexDouble );

    TCORE_zpotrf( clargs->uplo, clargs->n, tileA, &info );

    if ( (clargs->sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, clargs->sequence, clargs->request, clargs->iinfo+info );
    }
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zpotrf_cuda_func(void *descr[], void *cl_arg)
{
    cusolverDnHandle_t handle = starpu_cusolverDn_get_local_handle();
    struct cl_zpotrf_args_s *clargs = (struct cl_zpotrf_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;
    cuDoubleComplex *dW;
    int info, lwork;
    int *d_info;

    tileA  = cti_interface_get(descr[0]);
    tileW  = cti_interface_get(descr[1]);
    dW     = tileW->mat;
    lwork  = (tileW->ld - sizeof(int) ) / sizeof(cuDoubleComplex);
    d_info = (int *)(dW + lwork);

    TCUDA_zpotrf( clargs->uplo, clargs->n, tileA, dW, lwork, d_info, handle );

    cudaStreamSynchronize( starpu_cuda_get_local_stream() );

    cudaMemcpy( &info, d_info, sizeof(int), cudaMemcpyDeviceToHost );

    if ( ( clargs->sequence->status == CHAMELEON_SUCCESS ) && ( info != 0 ) ) {
        RUNTIME_sequence_flush( NULL, clargs->sequence, clargs->request, clargs->iinfo + info );
    }
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP) && defined(STARPU_HAVE_LIBHIPSOLVER)
static void
cl_zpotrf_hip_func( void *descr[], void *cl_arg )
{
    hipsolverDnHandle_t handle = starpu_hipsolverDn_get_local_handle();
    struct cl_zpotrf_args_s *clargs = (struct cl_zpotrf_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;
    hipDoubleComplex *dW;
    int info, lwork;
    int *d_info;

    tileA  = cti_interface_get(descr[0]);
    tileW  = cti_interface_get(descr[1]);
    dW     = tileW->mat;
    lwork  = (tileW->ld - sizeof(int) ) / sizeof(hipDoubleComplex);
    d_info = (int *)(dW + lwork);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );

    HIP_zpotrf( clargs->uplo, clargs->n, tileA->mat, tileA->ld, dW, lwork, d_info, handle );

    hipStreamSynchronize(starpu_hip_get_local_stream());

    hipMemcpy( &info, d_info, sizeof(int), hipMemcpyDeviceToHost);

    if ( (clargs->sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, clargs->sequence, clargs->request, clargs->iinfo+info );
    }
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP) && defined(STARPU_HAVE_LIBHIPSOLVER)
CODELETS_GPU( zpotrf, cl_zpotrf_cpu_func, cl_zpotrf_hip_func,  )
#else
CODELETS( zpotrf, cl_zpotrf_cpu_func, cl_zpotrf_cuda_func,  )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo )
{
    struct cl_zpotrf_args_s *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "zpotrf";
    CHAM_tile_t             *tileA;
    int                      is_rectask = 0;
    rectask_args_t          *rtargs     = NULL;
    (void)rtargs;

    /*
     * Add the workspace as scratch data when gpus are used
     */
    /* int accessWS = ( starpu_cuda_worker_get_count() | starpu_hip_worker_get_count() ) */
    /*                              ? ( STARPU_SCRATCH | STARPU_NOFOOTPRINT ) */
    /*                              : ( STARPU_NONE ); */
    int accessWS = options->ws_worker
        ? ( STARPU_SCRATCH | STARPU_NOFOOTPRINT )
        : ( STARPU_NONE );

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    tileA = A->get_blktile( A, Am, An );

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( tileA->format & CHAMELEON_TILE_DESC );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority;
        rtargs->tiles[0] = tileA;
        cl_name = "zpotrf_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zpotrf_args_s ) );
        clargs->uplo     = uplo;
        clargs->n        = n;
        clargs->iinfo    = iinfo;
        clargs->sequence = options->sequence;
        clargs->request  = options->request;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1, tileA );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zpotrf,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zpotrf_args_s),

        /* Task handles */
        STARPU_RW, RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessWS,  options->ws_worker,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zpotrf ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zpotrf( n ),

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zpotrf )
        0 );

    (void)tileA;
    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo )
{
    INSERT_TASK_COMMON_PARAMETERS( zpotrf, 2 );

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, A->get_rankof( A, Am, An ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_RW );

    /*
     * Add the workspace as scratch data when gpus are used
     */
    if ( options->ws_worker ) {
        starpu_cham_register_descr( &nbdata, descrs, options->ws_worker, STARPU_SCRATCH | STARPU_NOFOOTPRINT );
    }

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
        clargs = malloc( sizeof( struct cl_zpotrf_args_s ) );
        clargs->uplo     = uplo;
        clargs->n        = n;
        clargs->iinfo    = iinfo;
        clargs->sequence = options->sequence;
        clargs->request  = options->request;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zpotrf_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zpotrf_callback );

        /* Flops */
        task->flops = flops_zpotrf( n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 1,
                                             A->get_blktile( A, Am, An ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zpotrf", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
