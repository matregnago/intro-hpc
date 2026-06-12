/**
 *
 * @file starpu/codelet_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Gwenole Lucas
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztradd_args_s {
    cham_uplo_t           uplo;
    cham_trans_t          trans;
    int                   m;
    int                   n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_ztradd_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_ztradd_args_s *clargs  = (struct cl_ztradd_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pztradd( clargs->uplo, clargs->trans,
                       clargs->alpha, rtargs->tiles[0]->mat,
                       clargs->beta,  rtargs->tiles[1]->mat,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztradd_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztradd_args_s *clargs = (struct cl_ztradd_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_ztradd( clargs->uplo, clargs->trans, clargs->m, clargs->n,
                  clargs->alpha, tileA, clargs->beta, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( ztradd, cl_ztradd_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    struct cl_ztradd_args_s     *clargs  = NULL;
    int                          exec    = 0;
    const char                  *cl_name = "ztradd";
    enum starpu_data_access_mode accessB = STARPU_RW;
    CHAM_tile_t                 *tileA;
    CHAM_tile_t                 *tileB;
    int                          is_rectask = 0;
    rectask_args_t              *rtargs     = NULL;
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
        cl_name = "ztradd_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_ztradd_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Reduce the B access if possible */
    if ( ( uplo == ChamUpperLower ) && ( beta == 0. ) ) {
        accessB = STARPU_W;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztradd,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztradd_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,       RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztradd ),
        STARPU_NAME,              cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( ztradd )
        0 );

    (void)tileA;
    (void)tileB;
    (void)nb;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( ztradd, 2 );
    enum starpu_data_access_mode accessB = STARPU_RW;

    /* Reduce the B access if possible */
    if ( ( uplo == ChamUpperLower ) && ( beta == 0. ) ) {
        accessB = STARPU_W;
    }

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
        clargs = malloc( sizeof( struct cl_ztradd_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_ztradd_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_ztradd_callback );

        /* Flops */
        //task->flops = flops_ztradd( m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_ztradd", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
