/**
 *
 * @file starpu/codelet_ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm StarPU codelet
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
 * @author Terry Cojean
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztrsm_args_s {
    cham_side_t           side;
    cham_uplo_t           uplo;
    cham_trans_t          transA;
    cham_diag_t           diag;
    int                   m;
    int                   n;
    CHAMELEON_Complex64_t alpha;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_ztrsm_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_ztrsm_args_s *clargs  = (struct cl_ztrsm_args_s *)(t->cl_arg);
    rectask_args_t         *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t       request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pztrsm( clargs->side, clargs->uplo, clargs->transA, clargs->diag,
                      clargs->alpha, rtargs->tiles[0]->mat, rtargs->tiles[1]->mat,
                      rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztrsm_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztrsm_args_s *clargs = (struct cl_ztrsm_args_s*)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( tileA->flttype == ChamComplexDouble );
    assert( tileB->flttype == ChamComplexDouble );

    TCORE_ztrsm( clargs->side, clargs->uplo, clargs->transA, clargs->diag,
                 clargs->m, clargs->n, clargs->alpha, tileA, tileB );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_ztrsm_cuda_func(void *descr[], void *cl_arg)
{
    struct cl_ztrsm_args_s *clargs = (struct cl_ztrsm_args_s*)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCUDA_ztrsm( clargs->side, clargs->uplo, clargs->transA, clargs->diag,
                 clargs->m, clargs->n, (cuDoubleComplex *)&(clargs->alpha), tileA, tileB, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_ztrsm_hip_func(void *descr[], void *cl_arg)
{
    struct cl_ztrsm_args_s *clargs = (struct cl_ztrsm_args_s*)cl_arg;
    hipblasHandle_t         handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    HIP_ztrsm(
        clargs->side, clargs->uplo, clargs->transA, clargs->diag,
        clargs->m, clargs->n,
        (hipDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */

#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( ztrsm, cl_ztrsm_cpu_func, cl_ztrsm_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( ztrsm, cl_ztrsm_cpu_func, cl_ztrsm_cuda_func, STARPU_CUDA_ASYNC )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_ztrsm_args_s  *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "ztrsm";
    CHAM_tile_t             *tileA;
    CHAM_tile_t             *tileB;
    int                      is_rectask = 0;
    rectask_args_t          *rtargs     = NULL;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    tileA = A->get_blktile( A, Am, An );
    tileB = B->get_blktile( B, Bm, Bn );

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( ( tileA->format & CHAMELEON_TILE_DESC ) &&
                   ( tileB->format & CHAMELEON_TILE_DESC ) );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + sizeof(CHAM_tile_t*) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority ;
        rtargs->tiles[0] = tileA;
        rtargs->tiles[1] = tileB;
        cl_name = "ztrsm_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_ztrsm_args_s ) );
        clargs->side   = side;
        clargs->uplo   = uplo;
        clargs->transA = transA;
        clargs->diag   = diag;
        clargs->m      = m;
        clargs->n      = n;
        clargs->alpha  = alpha;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2, tileA, tileB );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztrsm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztrsm_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,     RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztrsm ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_ztrsm( side, m, n ),

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( ztrsm )
        0 );

    (void)tileA;
    (void)tileB;
    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn )
{
    INSERT_TASK_COMMON_PARAMETERS( ztrsm, 2 );

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_R  );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, B, Bm, Bn, STARPU_RW );

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
        clargs = malloc( sizeof( struct cl_ztrsm_args_s ) );
        clargs->side   = side;
        clargs->uplo   = uplo;
        clargs->transA = transA;
        clargs->diag   = diag;
        clargs->m      = m;
        clargs->n      = n;
        clargs->alpha  = alpha;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_ztrsm_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_ztrsm_callback );

        /* Flops */
        task->flops = flops_ztrsm( side, m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_ztrsm", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
