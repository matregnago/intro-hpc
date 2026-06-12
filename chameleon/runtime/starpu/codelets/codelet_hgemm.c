/**
 *
 * @file starpu/codelet_hgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hgemm StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelets.h"

CHAMELEON_CL_CB( hgemm, cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), cti_handle_get_n(task->handles[0]), 2. *M*N*K) /* If A^t, computation is wrong */

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
static void
cl_hgemm_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_hgemm_args_s *clargs = (struct cl_hgemm_args_s *)cl_arg;
    cublasHandle_t          handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    TCUDA_hgemm( clargs->transA, clargs->transB,
                 clargs->m, clargs->n, clargs->k,
                 (CHAMELEON_Real16_t*)&(clargs->alpha), tileA, tileB,
                 (CHAMELEON_Real16_t*)&(clargs->beta),  tileC,
                 handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_hgemm_hip_func( void *descr[], void *cl_arg )
{
    struct cl_hgemm_args_s *clargs = (struct cl_hgemm_args_s *)cl_arg;
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

    HIP_hgemm(
        clargs->transA, clargs->transB,
        clargs->m, clargs->n, clargs->k,
        (CHAMELEON_Real16_t*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        (CHAMELEON_Real16_t*)&(clargs->beta),
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
CODELETS_GPU( hgemm, NULL, cl_hgemm_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( hgemm, NULL, cl_hgemm_cuda_func, STARPU_CUDA_ASYNC )
#endif

void INSERT_TASK_hgemm_Astat( const RUNTIME_option_t *options,
                              cham_trans_t transA, cham_trans_t transB,
                              int m, int n, int k, int nb,
                              CHAMELEON_Real16_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                        const CHAM_desc_t *B, int Bm, int Bn,
                              CHAMELEON_Real16_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* if ( alpha == 0. ) { */
    /*     INSERT_TASK_hlascal( options, ChamUpperLower, m, n, nb, */
    /*                          beta, C, Cm, Cn ); */
    /*     return; */
    /* } */

    struct cl_hgemm_args_s  *clargs = NULL;
    int                      accessC;
    int                      exec    = 0;
    const char              *cl_name = "hgemm_Astat";

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
        clargs = malloc( sizeof( struct cl_hgemm_args_s ) );
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

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_hgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_hgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamRealHalf, Am, An),
        STARPU_R, RTBLKADDR(B, ChamRealHalf, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamRealHalf, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( hgemm ),
        STARPU_EXECUTE_ON_NODE,   A->get_rankof(A, Am, An),
        STARPU_NAME,              cl_name,
        0 );
}

void INSERT_TASK_hgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Real16_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                  const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Real16_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    /* if ( alpha == 0. ) { */
    /*     INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb, */
    /*                          beta, C, Cm, Cn ); */
    /*     return; */
    /* } */

    struct cl_hgemm_args_s  *clargs = NULL;
    int                      accessC;
    int                      exec = 0;
    const char              *cl_name = "hgemm";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_hgemm_args_s ) );
        clargs->transA = transA;
        clargs->transB = transB;
        clargs->m      = m;
        clargs->n      = n;
        clargs->k      = k;
        clargs->alpha  = alpha;
        clargs->beta   = beta;
    }

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : (STARPU_RW | ((beta == 1.) ? STARPU_COMMUTE : 0));

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_hgemm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_hgemm_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamRealHalf, Am, An),
        STARPU_R, RTBLKADDR(B, ChamRealHalf, Bm, Bn),
        accessC,  RTBLKADDR(C, ChamRealHalf, Cm, Cn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( hgemm ),
        STARPU_NAME,              cl_name,
        0 );
}
