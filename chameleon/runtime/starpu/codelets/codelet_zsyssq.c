/**
 *
 * @file starpu/codelet_zsyssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyssq StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zsyssq_cpu_func(void *descr[], void *cl_arg)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    tileA = cti_interface_get(descr[0]);
    tileW = cti_interface_get(descr[1]);
    starpu_codelet_unpack_args( cl_arg, &storev, &uplo, &n );
    TCORE_zsyssq( storev, uplo, n, tileA, tileW );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zsyssq, cl_zsyssq_cpu_func)

void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(SCALESUMSQ, SCALESUMSQm, SCALESUMSQn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zsyssq,
        STARPU_VALUE,    &storev,                     sizeof(cham_store_t),
        STARPU_VALUE,    &uplo,                       sizeof(int),
        STARPU_VALUE,    &n,                          sizeof(int),
        STARPU_R,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,       RTBLKADDR(SCALESUMSQ, ChamRealDouble, SCALESUMSQm, SCALESUMSQn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgessq ),
        0 );
}
