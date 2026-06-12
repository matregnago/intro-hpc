/**
 *
 * @file starpu/codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Gwenole Lucas
 * @author Brieuc Nicolas
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaset_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zlaset_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zlaset_args_s *clargs  = (struct cl_zlaset_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzlaset( clargs->uplo, clargs->alpha, clargs->beta, rtargs->tiles[0]->mat,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaset_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaset_args_s *clargs = (struct cl_zlaset_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlaset( clargs->uplo, clargs->m, clargs->n, clargs->alpha, clargs->beta, tileA );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zlaset_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_zlaset_args_s *clargs = (struct cl_zlaset_args_s *)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCUDA_zlaset( clargs->uplo, clargs->m, clargs->n,
                  (const cuDoubleComplex *)&(clargs->alpha),
                  (const cuDoubleComplex *)&(clargs->beta),
                  tileA, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zlaset, cl_zlaset_cpu_func, cl_zlaset_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *A, int Am, int An )
{
    struct cl_zlaset_args_s     *clargs  = NULL;
    int                          exec    = 0;
    const char                  *cl_name = "zlaset";
    enum starpu_data_access_mode accessA = STARPU_W;
    CHAM_tile_t                 *tileA;
    int                          is_rectask = 0;
    rectask_args_t              *rtargs     = NULL;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
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
        cl_name = "zlaset_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zlaset_args_s ) );
        clargs->uplo  = uplo;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    if ( uplo != ChamUpperLower ) {
        accessA = STARPU_RW;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlaset,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlaset_args_s),
        accessA, RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zlaset ),
        STARPU_NAME,              cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zlaset )
        0 );

    (void)tileA;
}
