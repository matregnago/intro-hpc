/**
 *
 * @file starpu/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm StarPU codelet
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
 * @author Philippe Swartvagher
 * @author Lucas Nesi
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

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zgemm_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zgemm_args_s *clargs  = (struct cl_zgemm_args_s *)(t->cl_arg);
    rectask_args_t         *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t       request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzgemm( NULL, clargs->transA, clargs->transB,
                      clargs->alpha, rtargs->tiles[0]->mat, rtargs->tiles[1]->mat,
                      clargs->beta,  rtargs->tiles[2]->mat,
                      rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgemm_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    assert( tileA->flttype == ChamComplexDouble );
    assert( tileB->flttype == ChamComplexDouble );
    assert( tileC->flttype == ChamComplexDouble );

    TCORE_zgemm( clargs->transA, clargs->transB,
                 clargs->m, clargs->n, clargs->k,
                 clargs->alpha, tileA, tileB,
                 clargs->beta,  tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zgemm_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    TCUDA_zgemm( clargs->transA, clargs->transB,
                 clargs->m, clargs->n, clargs->k,
                 (cuDoubleComplex *)&(clargs->alpha), tileA, tileB,
                 (cuDoubleComplex *)&(clargs->beta),  tileC, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_zgemm_hip_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s *clargs = (struct cl_zgemm_args_s *)cl_arg;
    hipblasHandle_t         handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );

    HIP_zgemm(
        clargs->transA, clargs->transB,
        clargs->m, clargs->n, clargs->k,
        (hipDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        (hipDoubleComplex*)&(clargs->beta),
        tileC->mat, tileC->ld,
        handle );

    return;
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zgemm, cl_zgemm_cpu_func, cl_zgemm_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zgemm, cl_zgemm_cpu_func, cl_zgemm_cuda_func, STARPU_CUDA_ASYNC )
#endif

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                              cham_trans_t transA, cham_trans_t transB,
                              int m, int n, int k, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                           const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zgemm_args_s  *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "zgemm_Astat";
    uint32_t                 where   = cl_zgemm.where;
    int                      accessC;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
     /* Check A as write, since it will be the owner of the computation */
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
     /* Check C as read, since it will be used in a reduction */
    CHAMELEON_ACCESS_R(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;
    }

    /* Reduce the C access if needed */
    if ( beta == 0. ) {
        accessC = STARPU_W;
    }
#if defined(HAVE_STARPU_MPI_REDUX)
    else if ( beta == 1. ) {
        accessC = STARPU_MPI_REDUX;
    }
#endif
    else {
        accessC = STARPU_RW;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* WARNING: CUDA 12.3 has an issue when m or n or k=1 in double complex,
       thus we disable gemm on gpu in these cases */
#if defined(PRECISION_z)
    if ( (k == 1) || (n == 1) || (m == 1) ) {
        where = STARPU_CPU;
    }
#endif

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R, RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgemm ),
        STARPU_EXECUTE_ON_NODE,   A->get_rankof(A, Am, An),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zgemm( m, n, k ),
        STARPU_EXECUTE_WHERE,     where,
        0 );
}

void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == (CHAMELEON_Complex64_t)0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zgemm_args_s  *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "zgemm";
    uint32_t                 where   = cl_zgemm.where;
    int                      accessC;
    CHAM_tile_t             *tileA;
    CHAM_tile_t             *tileB;
    CHAM_tile_t             *tileC;
    int                      is_rectask = 0;
    rectask_args_t          *rtargs     = NULL;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    tileA = A->get_blktile( A, Am, An );
    tileB = B->get_blktile( B, Bm, Bn );
    tileC = C->get_blktile( C, Cm, Cn );

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( ( tileA->format & CHAMELEON_TILE_DESC ) &&
                   ( tileB->format & CHAMELEON_TILE_DESC ) &&
                   ( tileC->format & CHAMELEON_TILE_DESC ) );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + 2 * sizeof(CHAM_tile_t*) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority;
        rtargs->tiles[0] = tileA;
        rtargs->tiles[1] = tileB;
        rtargs->tiles[2] = tileC;
        cl_name = "zgemm_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;
    }

    /* Reduce the C access if needed */
    if ( beta == (CHAMELEON_Complex64_t)0. ) {
        accessC = STARPU_W;
    }
    else {
        accessC = STARPU_RW;
        /* rectask is not compatible with commute yet */
        if ( (beta == (CHAMELEON_Complex64_t)1.) && !is_rectask ) {
            accessC |= STARPU_COMMUTE;
        }
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3, tileA, tileB, tileC );

    /* WARNING: CUDA 12.3 has an issue when m or n or k=1 in double complex,
       thus we disable gemm on gpu in these cases */
#if defined(PRECISION_z)
    if ( (k == 1) || (n == 1) || (m == 1) ) {
        where = STARPU_CPU;
    }
#endif

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R, RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgemm ),
        STARPU_NAME,              cl_name,
        STARPU_FLOPS,             flops_zgemm( m, n, k ),
        STARPU_EXECUTE_WHERE,     where,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zgemm )
        0 );

    (void)tileA;
    (void)tileB;
    (void)tileC;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void __INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                          int xrank, int accessC,
                          cham_trans_t transA, cham_trans_t transB,
                          int m, int n, int k, int nb,
                          CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                       const CHAM_desc_t *B, int Bm, int Bn,
                          CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == (CHAMELEON_Complex64_t)0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zgemm, 3 );

    /*
     * Register the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, xrank );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, A, Am, An, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, B, Bm, Bn, STARPU_R );
    starpu_cham_exchange_tile_before_execution( options, &params, &nbdata, descrs, C, Cm, Cn, accessC  );

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
#if defined(CHAMELEON_USE_MPI)
        assert( params.exchange_needed == 0 );
#endif
        return;
    }

    if ( params.do_execute )
    {
        int ret;
        struct starpu_task *task = starpu_task_create();
        task->cl = cl;

    /* WARNING: CUDA 12.3 has an issue when m or n or k=1 in double complex,
       thus we disable gemm on gpu in these cases */
#if defined(PRECISION_z)
        if ( (k == 1) || (n == 1) || (m == 1) ) {
            task->where = STARPU_CPU;
        }
#endif

        /* Set codelet parameters */
        clargs = malloc( sizeof( struct cl_zgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zgemm_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgemm_callback );

        /* Flops */
        task->flops = flops_zgemm( m, n, k );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 3,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ),
                                             C->get_blktile( C, Cm, Cn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zgemm", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

void INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                              cham_trans_t transA, cham_trans_t transB,
                              int m, int n, int k, int nb,
                              CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                           const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* Reduce the C access if needed */
    int accessC = ( beta == (CHAMELEON_Complex64_t)0. ) ? STARPU_W : STARPU_RW;

#if defined(HAVE_STARPU_MPI_REDUX)
    if ( beta == (CHAMELEON_Complex64_t)1. ) {
        accessC = STARPU_MPI_REDUX;
    }
#endif

    __INSERT_TASK_zgemm( options,
                         A->get_rankof( A, Am, An ), accessC,
                         transA, transB, m, n, k, nb,
                         alpha, A, Am, An,
                                B, Bm, Bn,
                         beta,  C, Cm, Cn );
}

void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* Reduce the C access if needed */
    int accessC = ( beta == (CHAMELEON_Complex64_t)0. ) ? STARPU_W :
        (STARPU_RW | ((beta == (CHAMELEON_Complex64_t)1.) ? STARPU_COMMUTE : 0));

    __INSERT_TASK_zgemm( options,
                         C->get_rankof( C, Cm, Cn ), accessC,
                         transA, transB, m, n, k, nb,
                         alpha, A, Am, An,
                                B, Bm, Bn,
                         beta,  C, Cm, Cn );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
