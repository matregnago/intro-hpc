/**
 *
 * @file quark/codelet_ztrasm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrasm Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_ztrasm_quark(Quark *quark)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    cham_diag_t diag;
    int M;
    int N;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    quark_unpack_args_7( quark, storev, uplo, diag, M, N, tileA, tileW );
    TCORE_ztrasm( storev, uplo, diag, M, N, tileA, tileW->mat );
}

void INSERT_TASK_ztrasm(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_ztrasm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &storev,    VALUE,
        sizeof(int),              &uplo,      VALUE,
        sizeof(int),              &diag,      VALUE,
        sizeof(int),                     &M,         VALUE,
        sizeof(int),                     &N,         VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),  INPUT,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn),                 INOUT,
        0);
}
