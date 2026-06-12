/**
 *
 * @file quark/codelet_zunmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
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

void CORE_zunmlq_quark(Quark *quark)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_11(quark, side, trans, m, n, k, ib,
                         tileA, tileT, tileC, WORK, ldwork);
    TCORE_zunmlq(side, trans, m, n, k, ib,
                tileA, tileT, tileC, WORK, ldwork);
}

void INSERT_TASK_zunmlq(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m, int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_UNMLQ;
    QUARK_Insert_Task(opt->quark, CORE_zunmlq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(int),                     &k,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT | QUARK_REGION_U,
        sizeof(void*), RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), INPUT,
        sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn), INOUT,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL,      SCRATCH,
        sizeof(int),                     &nb,    VALUE,
        0);
}
