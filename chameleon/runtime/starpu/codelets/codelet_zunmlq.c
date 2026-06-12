/**
 *
 * @file starpu/codelet_zunmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zunmlq_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m, &n, &k, &ib, &ldW );

    TCORE_zunmlq( side, trans, m, n, k, ib,
                  tileA, tileT, tileC, tileW->mat, ldW );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zunmlq_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m, &n, &k, &ib, &ldW );

    TCUDA_zunmlqt( side, trans, m, n, k, ib,
                   tileA, tileT, tileC, tileW, ldW, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zunmlq, cl_zunmlq_cpu_func, cl_zunmlq_cuda_func, STARPU_CUDA_ASYNC)

void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn )
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zunmlq,
        STARPU_VALUE,    &side,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,         RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_RW,        RTBLKADDR(C, ChamComplexDouble, Cm, Cn),
        /* ib * nb */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &nb,                sizeof(int),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zunmlq ),
        0 );
}
