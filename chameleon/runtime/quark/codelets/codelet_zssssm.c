/**
 *
 * @file quark/codelet_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm Quark codelet
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

void CORE_zssssm_quark(Quark *quark)
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

    quark_unpack_args_11(quark, m1, n1, m2, n2, k, ib, tileA1, tileA2, tileL1, tileL2, IPIV);
    TCORE_zssssm(m1, n1, m2, n2, k, ib, tileA1, tileA2, tileL1, tileL2, IPIV);
}

void INSERT_TASK_zssssm(const RUNTIME_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n,
                       const CHAM_desc_t *A2, int A2m, int A2n,
                       const CHAM_desc_t *L1, int L1m, int L1n,
                       const CHAM_desc_t *L2, int L2m, int L2n,
                       const int *IPIV)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_SSSSM;
    QUARK_Insert_Task(opt->quark, CORE_zssssm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(void*), RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n),            INOUT,
        sizeof(void*), RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n),            INOUT | LOCALITY,
        sizeof(void*), RTBLKADDR(L1, CHAMELEON_Complex64_t, L1m, L1n),            INPUT,
        sizeof(void*), RTBLKADDR(L2, CHAMELEON_Complex64_t, L2m, L2n),            INPUT,
        sizeof(int)*nb,                      IPIV,          INPUT,
        0);
}
