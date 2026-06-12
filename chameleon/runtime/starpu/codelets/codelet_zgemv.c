/**
 *
 * @file starpu/codelet_zgemv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemv StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgemv_cpu_func(void *descr[], void *cl_arg)
{
    cham_trans_t trans;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileX;
    int incX;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileY;
    int incY;

    tileA = cti_interface_get(descr[0]);
    tileX = cti_interface_get(descr[1]);
    tileY = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &trans, &m, &n, &alpha, &incX, &beta, &incY );
    TCORE_zgemv( trans, m, n,
                 alpha, tileA, tileX, incX,
                 beta,  tileY, incY );
}

#if defined(CHAMELEON_USE_CUDA) & 0
static void cl_zgemv_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    cuDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    cuDoubleComplex beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &beta);

    CUDA_zgemv(
        transA, transB,
        m, n, k,
        &alpha, tileA->mat, tileA->ld,
                tileB->mat, tileB->ld,
        &beta,  tileC->mat, tileC->ld,
        handle );
    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgemv, cl_zgemv_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY )
{
    int accessY = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(X, Xm, Xn);
    CHAMELEON_ACCESS_RW(Y, Ym, Yn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgemv,
        STARPU_VALUE,    &trans,             sizeof(cham_trans_t),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &alpha,             sizeof(CHAMELEON_Complex64_t),
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,         RTBLKADDR(X, ChamComplexDouble, Xm, Xn),
        STARPU_VALUE,    &incX,              sizeof(int),
        STARPU_VALUE,    &beta,              sizeof(CHAMELEON_Complex64_t),
        accessY,          RTBLKADDR(Y, ChamComplexDouble, Ym, Yn),
        STARPU_VALUE,    &incY,              sizeof(int),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgemv ),
        0);
}
