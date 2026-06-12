/**
 *
 * @file starpu/codelet_dzasum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dzasum StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
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
static void cl_dzasum_cpu_func(void *descr[], void *cl_arg)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tilework;

    tileA    = cti_interface_get(descr[0]);
    tilework = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &storev, &uplo, &M, &N);
    TCORE_dzasum(storev, uplo, M, N, tileA, tilework->mat );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dzasum, cl_dzasum_cpu_func)

void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_dzasum,

        /* Task codelet arguments */
        STARPU_VALUE, &storev, sizeof(cham_store_t),
        STARPU_VALUE, &uplo,   sizeof(cham_uplo_t),
        STARPU_VALUE, &M,      sizeof(int),
        STARPU_VALUE, &N,      sizeof(int),

        /* Task handles */
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,        RTBLKADDR(B, ChamRealDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( dzasum ),
        0);
}
