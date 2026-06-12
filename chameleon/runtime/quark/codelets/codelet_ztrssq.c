/**
 *
 * @file quark/codelet_ztrssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrssq Quark codelet
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

void CORE_ztrssq_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_diag_t diag;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    quark_unpack_args_6( quark, uplo, diag, m, n, tileA, tileW );
    TCORE_ztrssq( uplo, diag, m, n, tileA, tileW );
}

void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_diag_t diag,
                        int m, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TRSSQ;
    QUARK_Insert_Task(opt->quark, CORE_ztrssq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo, VALUE,
        sizeof(int),              &diag, VALUE,
        sizeof(int),                     &m,    VALUE,
        sizeof(int),                     &n,    VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(void*), RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INOUT,
        0);
}
