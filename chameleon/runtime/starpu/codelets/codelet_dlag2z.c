/**
 *
 * @file starpu/codelet_dlag2z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlag2z StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-10-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_dlag2z_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &m, &n);
    TCORE_dlag2z( uplo, m, n, tileA, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dlag2z, cl_dlag2z_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    enum starpu_data_access_mode accessB = STARPU_W;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( uplo != ChamUpperLower ) {
        accessB = STARPU_RW;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_dlag2z,

        /* Task codelet arguments */
        STARPU_VALUE, &uplo, sizeof(cham_uplo_t),
        STARPU_VALUE, &m,    sizeof(int),
        STARPU_VALUE, &n,    sizeof(int),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamRealDouble, Am, An),
        accessB,  RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( dlag2z ),
        0);
}
