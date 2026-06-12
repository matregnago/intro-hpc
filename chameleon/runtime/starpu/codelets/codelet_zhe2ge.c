/**
 *
 * @file starpu/codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge StarPU codelet
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
static void cl_zhe2ge_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N);
    TCORE_zhe2ge(uplo, M, N, tileA, tileB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhe2ge, cl_zhe2ge_cpu_func)

void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    (void)mb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zhe2ge,
        STARPU_VALUE,  &uplo,                sizeof(int),
        STARPU_VALUE,     &m,                        sizeof(int),
        STARPU_VALUE,     &n,                        sizeof(int),
        STARPU_R,             RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,             RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zhe2ge ),
        0 );
}
