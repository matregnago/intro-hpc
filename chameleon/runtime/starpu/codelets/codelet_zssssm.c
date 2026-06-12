/**
 *
 * @file starpu/codelet_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm StarPU codelet
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
static void cl_zssssm_cpu_func(void *descr[], void *cl_arg)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAM_tile_t *tileA1;
    CHAM_tile_t *tileA2;
    CHAM_tile_t *tileL1;
    CHAM_tile_t *tileL2;
    int *IPIV;

    tileA1 = cti_interface_get(descr[0]);
    tileA2 = cti_interface_get(descr[1]);
    tileL1 = cti_interface_get(descr[2]);
    tileL2 = cti_interface_get(descr[3]);
    starpu_codelet_unpack_args(cl_arg, &m1, &n1, &m2, &n2, &k, &ib, &IPIV);
    TCORE_zssssm(m1, n1, m2, n2, k, ib, tileA1, tileA2, tileL1, tileL2, IPIV);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zssssm, cl_zssssm_cpu_func)

void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV )
{
    (void)nb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A1, A1m, A1n);
    CHAMELEON_ACCESS_RW(A2, A2m, A2n);
    CHAMELEON_ACCESS_R(L1, L1m, L1n);
    CHAMELEON_ACCESS_R(L2, L2m, L2n);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zssssm,
        STARPU_VALUE,    &m1,                        sizeof(int),
        STARPU_VALUE,    &n1,                        sizeof(int),
        STARPU_VALUE,    &m2,                        sizeof(int),
        STARPU_VALUE,    &n2,                        sizeof(int),
        STARPU_VALUE,     &k,                        sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_RW,            RTBLKADDR(A1, ChamComplexDouble, A1m, A1n),
        STARPU_RW,            RTBLKADDR(A2, ChamComplexDouble, A2m, A2n),
        STARPU_R,             RTBLKADDR(L1, ChamComplexDouble, L1m, L1n),
        STARPU_R,             RTBLKADDR(L2, ChamComplexDouble, L2m, L2n),
        STARPU_VALUE,          &IPIV,                      sizeof(int*),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zssssm ),
        0 );
}
