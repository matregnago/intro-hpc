/**
 *
 * @file quark/codelet_zsyssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyssq Quark codelet
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

void CORE_zsyssq_quark(Quark *quark)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    quark_unpack_args_5( quark, storev, uplo, n, tileA, tileW );
    TCORE_zsyssq( storev, uplo, n, tileA, tileW );
}

void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                        cham_store_t storev, cham_uplo_t uplo, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_SYSSQ;
    QUARK_Insert_Task(opt->quark, CORE_zsyssq_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_store_t),            &storev, VALUE,
        sizeof(int),                     &uplo, VALUE,
        sizeof(int),                     &n,    VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(void*), RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INOUT,
        0);
}
