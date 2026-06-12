/**
 *
 * @file starpu/codelet_ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Lucas Barros de Assis
 * @author Samuel Thibault
 * @date 2024-10-21
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztpqrt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int L;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileWORK;

    tileA    = cti_interface_get(descr[0]);
    tileB    = cti_interface_get(descr[1]);
    tileT    = cti_interface_get(descr[2]);
    tileWORK = cti_interface_get(descr[3]); /* ib * nb */
    starpu_codelet_unpack_args( cl_arg, &M, &N, &L, &ib );

    TCORE_zlaset( ChamUpperLower, ib, N, 0., 0., tileT );
    TCORE_ztpqrt( M, N, L, ib,
                  tileA, tileB, tileT, tileWORK->mat );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztpqrt, cl_ztpqrt_cpu_func)

void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn )
{
    const char            *cl_name;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_ACCESS_W(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = (L == 0) ? "ztsqrt" : "zttqrt";

    rt_starpu_insert_task(
        &cl_ztpqrt,
        /* Task codelet arguments */
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,    sizeof(int),

        /* Task handles */
        STARPU_RW,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,      RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_W,       RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_SCRATCH, options->ws_worker,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztpqrt ),
        STARPU_EXECUTE_ON_NODE,   B->get_rankof(B, Bm, Bn),
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
