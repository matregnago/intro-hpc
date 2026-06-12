/**
 *
 * @file quark/codelet_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum Quark codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
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

void CORE_zlauum_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int N;
    CHAM_tile_t *tileA;

    quark_unpack_args_3(quark, uplo, N, tileA);
    TCORE_zlauum(uplo, N, tileA);
}

void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(opt->quark, CORE_zlauum_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,  VALUE,
        sizeof(int),                &n,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INOUT,
        0);

    (void)nb;
}
