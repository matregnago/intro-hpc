/**
 *
 * @file quark/codelet_zgelqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqt Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgelqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAMELEON_Complex64_t *TAU;
    CHAMELEON_Complex64_t *WORK;

    quark_unpack_args_7(quark, m, n, ib, tileA, tileT, TAU, WORK);
    TCORE_zlaset( ChamUpperLower, ib, m, 0., 0., tileT );
    TCORE_zgelqt(m, n, ib, tileA,  tileT,  TAU, WORK);
}

void INSERT_TASK_zgelqt(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GELQT;
    QUARK_Insert_Task(opt->quark, CORE_zgelqt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INOUT,
        sizeof(void*), RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), OUTPUT,
        sizeof(CHAMELEON_Complex64_t)*nb,    NULL,          SCRATCH,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL,          SCRATCH,
        0);
}
