/**
 *
 * @file starpu/codelet_ztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri StarPU codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Gwenole Lucas
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztrtri_args_s {
    cham_uplo_t         uplo;
    cham_diag_t         diag;
    int                 n;
    int                 iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_ztrtri_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_ztrtri_args_s *clargs  = (struct cl_ztrtri_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pztrtri( clargs->uplo, clargs->diag, rtargs->tiles[0]->mat,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztrtri_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztrtri_args_s *clargs = (struct cl_ztrtri_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    TCORE_ztrtri( clargs->uplo, clargs->diag, clargs->n, tileA, &info );

    if ( (clargs->sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, clargs->sequence, clargs->request, clargs->iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( ztrtri, cl_ztrtri_cpu_func )

void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo )
{
    struct cl_ztrtri_args_s *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "ztrtri";
    CHAM_tile_t             *tileA;
    int                      is_rectask = 0;
    rectask_args_t          *rtargs     = NULL;
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
        rtargs->tiles[0] = tileA;
        cl_name = "ztrtri_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_ztrtri_args_s ) );
        clargs->uplo     = uplo;
        clargs->diag     = diag;
        clargs->n        = n;
        clargs->iinfo    = iinfo;
        clargs->sequence = options->sequence;
        clargs->request  = options->request;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1, tileA );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztrtri,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztrtri_args_s),
        STARPU_RW,     RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztrtri ),
        STARPU_NAME,              cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( ztrtri )
        0 );

    (void)tileA;
    (void)nb;
}
