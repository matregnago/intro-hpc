/**
 *
 * @file starpu/codelet_zplgsy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgsy StarPU codelet
 *
 * @version 1.4.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
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

struct cl_zplgsy_args_s {
    CHAMELEON_Complex64_t  bump;
    cham_uplo_t            uplo;
    int                    m;
    int                    n;
    int                    bigM;
    int                    m0;
    int                    n0;
    unsigned long long int seed;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zplgsy_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zplgsy_args_s *clargs  = (struct cl_zplgsy_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzplgsy( clargs->bump, ChamUpperLower, rtargs->tiles[0]->mat,
                       clargs->bigM, clargs->m0, clargs->n0, clargs->seed,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplgsy_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zplgsy_args_s *clargs = (struct cl_zplgsy_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zplgsy( clargs->bump, clargs->m, clargs->n, tileA,
                  clargs->bigM, clargs->m0, clargs->n0, clargs->seed );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zplgsy, cl_zplgsy_cpu_func )

void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    struct cl_zplgsy_args_s *clargs = NULL;
    int                      exec = 0;
    const char              *cl_name = "zplgsy";
    CHAM_tile_t             *tileA;
    int                      is_rectask = 0;
    rectask_args_t          *rtargs     = NULL;
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
        rtargs->priority = options->priority;
        rtargs->tiles[0] = tileA;
        cl_name = "zplgsy_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zplgsy_args_s ) );
        clargs->bump = bump;
        clargs->uplo = ChamUpperLower;
        clargs->m    = m;
        clargs->n    = n;
        clargs->bigM = bigM;
        clargs->m0   = m0;
        clargs->n0   = n0;
        clargs->seed = seed;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zplgsy,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zplgsy_args_s),
        STARPU_W,      RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zplgsy ),
        STARPU_NAME,              cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zplgsy )
        0 );

    (void)tileA;
}
