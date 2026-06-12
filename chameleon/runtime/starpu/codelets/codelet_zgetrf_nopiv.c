/**
 *
 * @file starpu/codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv StarPU codelet
 *
 * @version 1.4.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgetrf_nopiv_args_s {
    int                 m;
    int                 n;
    int                 ib;
    int                 iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zgetrf_nopiv_rectask_func( struct starpu_task *t, void *_args )
{
    rectask_args_t   *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzgetrf_nopiv( NULL, rtargs->tiles[0]->mat, rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

/*
 * Codelet CPU
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_nopiv_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_nopiv_args_s *clargs = (struct cl_zgetrf_nopiv_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    TCORE_zgetrf_nopiv( clargs->m, clargs->n, clargs->ib, tileA, &info );

    if ( (clargs->sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, clargs->sequence, clargs->request, clargs->iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_nopiv, cl_zgetrf_nopiv_cpu_func)

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgetrf_nopiv(const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              int iinfo)
{
    struct cl_zgetrf_nopiv_args_s *clargs  = NULL;
    int                            exec    = 0;
    const char                    *cl_name = "zgetrf_nopiv";
    CHAM_tile_t                   *tileA;
    int                            is_rectask = 0;
    rectask_args_t                *rtargs     = NULL;
    (void)rtargs;

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
        rtargs->priority = options->priority ;
        rtargs->tiles[0] = tileA;
        cl_name = "zgetrf_nopiv_rectask";
    }
#endif

    /* Set codelet parameters */
    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zgetrf_nopiv_args_s ) );
        clargs->m        = m;
        clargs->n        = n;
        clargs->ib       = ib;
        clargs->iinfo    = iinfo;
        clargs->sequence = options->sequence;
        clargs->request  = options->request;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1, tileA );

    rt_starpu_insert_task(
        &cl_zgetrf_nopiv,
        /* Task codelet arguments */
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zgetrf_nopiv_args_s),
        STARPU_RW,                RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgetrf_nopiv ),
        STARPU_NAME,              cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zgetrf_nopiv )
        0 );

    (void)tileA;
    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_nopiv(const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              int iinfo)
{
    INSERT_TASK_COMMON_PARAMETERS( zgetrf_nopiv, 1 );

    /*
     * Register the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, A->get_rankof( A, Am, An ) );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_RW );

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
        clargs = malloc( sizeof( struct cl_zgetrf_nopiv_args_s ) );
        clargs->m        = m;
        clargs->n        = n;
        clargs->ib       = ib;
        clargs->iinfo    = iinfo;
        clargs->sequence = options->sequence;
        clargs->request  = options->request;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zgetrf_nopiv_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgetrf_nopiv_callback );

        /* Flops */
        task->flops = flops_zgetrf( m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 1,
                                             A->get_blktile( A, Am, An ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zgetrf_nopiv", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
