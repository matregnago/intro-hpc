/**
 *
 * @file starpu/codelet_zlanhe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlanhe_cpu_func(void *descr[], void *cl_arg)
{
    CHAM_tile_t *tilenormA;
    cham_normtype_t norm;
    cham_uplo_t uplo;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tilework;

    tileA     = cti_interface_get(descr[0]);
    tilework  = cti_interface_get(descr[1]);
    tilenormA = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &norm, &uplo, &N);
    TCORE_zlanhe( norm, uplo, N, tileA, tilework->mat, tilenormA->mat );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlanhe, cl_zlanhe_cpu_func)

void INSERT_TASK_zlanhe(const RUNTIME_option_t *options,
                       cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zlanhe,
        STARPU_VALUE,    &norm,              sizeof(int),
        STARPU_VALUE,    &uplo,              sizeof(int),
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_R,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_SCRATCH,  options->ws_worker,
        STARPU_W,        RTBLKADDR(B, ChamRealDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zlange ),
        0 );

    (void)NB;
}
