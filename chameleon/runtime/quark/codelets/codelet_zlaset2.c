/**
 *
 * @file quark/codelet_zlaset2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
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

void CORE_zlaset2_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;

    quark_unpack_args_5(quark, uplo, M, N, alpha, tileA);
    TCORE_zlaset2(uplo, M, N, alpha, tileA);
}

void INSERT_TASK_zlaset2(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int M, int N,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASET;
    QUARK_Insert_Task(opt->quark, CORE_zlaset2_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha, VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),      OUTPUT,
        0);
}
