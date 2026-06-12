/**
 *
 * @file starpu/codelet_zgessm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm StarPU codelet
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
static void cl_zgessm_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    CHAM_tile_t *tileD;
    CHAM_tile_t *tileA;

    tileD = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &k, &ib, &IPIV);
    TCORE_zgessm(m, n, k, ib, IPIV, tileD, tileA);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgessm, cl_zgessm_cpu_func)

void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An )
{
    (void)nb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(L, Lm, Ln);
    CHAMELEON_ACCESS_R(D, Dm, Dn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgessm,
        STARPU_VALUE,     &m,                        sizeof(int),
        STARPU_VALUE,     &n,                        sizeof(int),
        STARPU_VALUE,     &k,                        sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_VALUE,          &IPIV,                      sizeof(int*),
        STARPU_R,             RTBLKADDR(L, ChamComplexDouble, Lm, Ln),
        STARPU_R,             RTBLKADDR(D, ChamComplexDouble, Dm, Dn),
        STARPU_RW,            RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgessm ),
        0 );
}
