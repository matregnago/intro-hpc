/**
 *
 * @file quark/codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro Quark codelet
 *
 * @version 1.4.0
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zlatro_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int M;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    quark_unpack_args_6(quark, uplo, trans, M, N, tileA, tileB);
    TCORE_zlatro(uplo, trans, M, N, tileA, tileB);
}

void INSERT_TASK_zlatro(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zlatro_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), OUTPUT,
        0);

    (void)nb;
}
