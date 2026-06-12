/**
 *
 * @file quark/codelet_ztsmlq_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq_hetra1 Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_ztsmlq_hetra1_quark(Quark *quark)
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
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m1, n1, m2, n2, k, ib, tileA1, tileA2, tileV, tileT, WORK, ldwork);
    TCORE_ztsmlq_hetra1( side, trans, m1, n1, m2, n2, k, ib,
                         tileA1, tileA2, tileV, tileT, WORK, ldwork );
}

void INSERT_TASK_ztsmlq_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n,
                              const CHAM_desc_t *A2, int A2m, int A2n,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int ldwork = side == ChamLeft ? ib : nb;

    QUARK_Insert_Task(opt->quark, CORE_ztsmlq_hetra1_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,   VALUE,
        sizeof(int),              &trans,  VALUE,
        sizeof(int),                     &m1,     VALUE,
        sizeof(int),                     &n1,     VALUE,
        sizeof(int),                     &m2,     VALUE,
        sizeof(int),                     &n2,     VALUE,
        sizeof(int),                     &k,      VALUE,
        sizeof(int),                     &ib,     VALUE,
        sizeof(void*), RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n), INOUT|QUARK_REGION_U|QUARK_REGION_D,
        sizeof(void*), RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n), INOUT,
        sizeof(void*), RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn),    INPUT,
        sizeof(void*), RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),    INPUT,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL,   SCRATCH,
        sizeof(int),                     &ldwork, VALUE,
        0);
}
