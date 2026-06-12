/**
 *
 * @file starpu/codelet_zgeqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrt StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgeqrt_cpu_func(void *descr[], void *cl_arg)
{
    CHAMELEON_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileW;
    CHAMELEON_Complex64_t *TAU;
    CHAMELEON_Complex64_t *WORK;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileW = cti_interface_get(descr[2]); /* max(m,n) + ib * n */

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &h_work);

    TAU  = tileW->mat;
    WORK = TAU + chameleon_max( m, n );

    TCORE_zlaset( ChamUpperLower, ib, n, 0., 0., tileT );
    TCORE_zgeqrt(m, n, ib, tileA, tileT, TAU, WORK );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgeqrt, cl_zgeqrt_cpu_func)

void INSERT_TASK_zgeqrt(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn)
{
    (void)nb;
    CHAMELEON_starpu_ws_t *h_work = (CHAMELEON_starpu_ws_t*)(options->ws_host);

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_W(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgeqrt,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,         RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        /* max( nb * (ib+1), ib * (ib+nb) ) */
        STARPU_SCRATCH,   options->ws_worker,
        /* ib * (m+3*ib) + max(m,n) */
        STARPU_VALUE,    &h_work,            sizeof(CHAMELEON_starpu_ws_t *),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgeqrt ),
        0 );
}
