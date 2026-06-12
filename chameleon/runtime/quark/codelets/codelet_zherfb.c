/**
 *
 * @file quark/codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zherfb_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int n;
    int k;
    int ib;
    int nb;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_10(quark, uplo, n, k, ib, nb, tileA, tileT, tileC, WORK, ldwork);
    TCORE_zherfb(uplo, n, k, ib, nb, tileA, tileT, tileC, WORK, ldwork);
}

void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zherfb_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo, VALUE,
        sizeof(int),                       &n,    VALUE,
        sizeof(int),                       &k,    VALUE,
        sizeof(int),                       &ib,   VALUE,
        sizeof(int),                       &nb,   VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), (uplo == ChamUpper) ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(void*), RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), INPUT,
        sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn), (uplo == ChamUpper) ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(CHAMELEON_Complex64_t)*2*nb*nb,  NULL, SCRATCH,
        sizeof(int),                       &nb,   VALUE,
        0);
}
