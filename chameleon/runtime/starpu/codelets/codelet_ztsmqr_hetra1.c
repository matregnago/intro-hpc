/**
 *
 * @file starpu/codelet_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmqr_hetra1 StarPU codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
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
static void cl_ztsmqr_hetra1_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAM_tile_t *tileA1;
    CHAM_tile_t *tileA2;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileW;
    int ldW;

    tileA1 = cti_interface_get(descr[0]);
    tileA2 = cti_interface_get(descr[1]);
    tileV  = cti_interface_get(descr[2]);
    tileT  = cti_interface_get(descr[3]);
    tileW  = cti_interface_get(descr[4]);

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &ldW );
    TCORE_ztsmqr_hetra1( side, trans, m1, n1, m2, n2, k, ib,
                         tileA1, tileA2, tileV, tileT, tileW->mat, ldW );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztsmqr_hetra1, cl_ztsmqr_hetra1_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_ztsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V,  int Vm,  int Vn,
                                const CHAM_desc_t *T,  int Tm,  int Tn )
{

    int ldWORK = side == ChamLeft ? ib : nb;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A1, A1m, A1n);
    CHAMELEON_ACCESS_RW(A2, A2m, A2n);
    CHAMELEON_ACCESS_R(V, Vm, Vn);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_ztsmqr_hetra1,
        STARPU_VALUE,    &side,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &m1,                sizeof(int),
        STARPU_VALUE,    &n1,                sizeof(int),
        STARPU_VALUE,    &m2,                sizeof(int),
        STARPU_VALUE,    &n2,                sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_VALUE,    &ldWORK,            sizeof(int),
        STARPU_RW,        RTBLKADDR(A1, ChamComplexDouble, A1m, A1n),
        STARPU_RW,        RTBLKADDR(A2, ChamComplexDouble, A2m, A2n),
        STARPU_R,         RTBLKADDR(V, ChamComplexDouble, Vm, Vn),
        STARPU_R,         RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_SCRATCH,   options->ws_worker,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztsmqr_hetra1 ),
        0 );
}
