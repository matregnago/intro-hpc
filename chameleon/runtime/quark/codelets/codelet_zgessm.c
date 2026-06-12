/**
 *
 * @file quark/codelet_zgessm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgessm Quark codelet
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
#include "coreblas/cblas.h"
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgessm_quark(Quark *quark)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    CHAM_tile_t *tileL;
    CHAM_tile_t *tileD;
    CHAM_tile_t *tileA;

    quark_unpack_args_8(quark, m, n, k, ib, IPIV, tileL, tileD, tileA);
    TCORE_zgessm(m, n, k, ib, IPIV, tileD, tileA);
}

void INSERT_TASK_zgessm(const RUNTIME_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       const CHAM_desc_t *L, int Lm, int Ln,
                       const CHAM_desc_t *D, int Dm, int Dn,
                       const CHAM_desc_t *A, int Am, int An)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GESSM;
    QUARK_Insert_Task(opt->quark, CORE_zgessm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        sizeof(void*), RTBLKADDR(L, CHAMELEON_Complex64_t, Lm, Ln),             INPUT | QUARK_REGION_L,
        sizeof(void*), RTBLKADDR(D, CHAMELEON_Complex64_t, Dm, Dn),             INPUT | QUARK_REGION_L,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INOUT,
        0);
}
