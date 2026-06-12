/**
 *
 * @file starpu/codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zherfb_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n, k, ib, nb;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &uplo, &n, &k, &ib, &nb, &ldW );

    TCORE_zherfb( uplo, n, k, ib, nb, tileA, tileT, tileC, tileW->mat, ldW );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zherfb_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_uplo_t uplo;
    int n, k, ib, nb;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &uplo, &n, &k, &ib, &nb, &ldW );

    TCUDA_zherfb( uplo, n, k, ib, nb, tileA, tileT, tileC, tileW, ldW, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zherfb, cl_zherfb_cpu_func, cl_zherfb_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zherfb,
        STARPU_VALUE,    &uplo,              sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_VALUE,    &nb,                sizeof(int), /* ldw */
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,         RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_RW,        RTBLKADDR(C, ChamComplexDouble, Cm, Cn),
        STARPU_SCRATCH,   options->ws_worker,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zherfb ),
        0 );
}
