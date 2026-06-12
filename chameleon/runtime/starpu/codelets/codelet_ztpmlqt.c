/**
 *
 * @file starpu/codelet_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmlqt StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztpmlqt_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* ib * nb */
    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    TCORE_ztpmlqt( side, trans, M, N, K, L, ib,
                   tileV, tileT, tileA, tileB, tileW->mat );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_ztpmlqt_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* 3*ib*nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    TCUDA_ztpmlqt( side, trans, M, N, K, L, ib,
                   tileV, tileT, tileA, tileB, tileW, lwork, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(ztpmlqt, cl_ztpmlqt_cpu_func, cl_ztpmlqt_cuda_func, STARPU_CUDA_ASYNC)

void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn )
{
    const char            *cl_name;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(V, Vm, Vn);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = (L == 0) ? "ztsmlqt" : "zttmlqt";

    rt_starpu_insert_task(
        &cl_ztpmlqt,
        /* Task codelet arguments */
        STARPU_VALUE, &side,  sizeof(int),
        STARPU_VALUE, &trans, sizeof(int),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &K,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,     sizeof(int),
        STARPU_VALUE, &(options->ws_wsize), sizeof(size_t),

        /* Task handles */
        STARPU_R,      RTBLKADDR(V, ChamComplexDouble, Vm, Vn),
        STARPU_R,      RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_RW,     RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,     RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_SCRATCH,   options->ws_worker,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztpmlqt ),
        STARPU_EXECUTE_ON_NODE,   B->get_rankof(B, Bm, Bn),
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
