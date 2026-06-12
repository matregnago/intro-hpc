/**
 *
 * @file starpu/codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @author Lucas Barros de Assis
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zsytrf_nopiv_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int n;
    CHAM_tile_t *tileA;
    int iinfo;

    tileA = cti_interface_get(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &iinfo);
    TCORE_zsytf2_nopiv(uplo, n, tileA);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zsytrf_nopiv, cl_zsytrf_nopiv_cpu_func)

void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo )
{
    (void)nb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zsytrf_nopiv,
        STARPU_VALUE,    &uplo,                      sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_RW,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        /* STARPU_SCRATCH,   options->ws_worker, */

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zsytrf_nopiv ),
        0 );
}
