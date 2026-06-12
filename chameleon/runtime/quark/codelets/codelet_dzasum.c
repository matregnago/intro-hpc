/**
 *
 * @file quark/codelet_dzasum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dzasum Quark codelet
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

static inline void
CORE_dzasum_quark(Quark *quark)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAM_tile_t *A;
    CHAM_tile_t *work;

    quark_unpack_args_6(quark, storev, uplo, M, N, A, work);
    TCORE_dzasum( storev, uplo, M, N, A, work->mat );
}

void INSERT_TASK_dzasum(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, int M, int N,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_ASUM;
    QUARK_Insert_Task(opt->quark, CORE_dzasum_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &storev,    VALUE,
        sizeof(int),              &uplo,      VALUE,
        sizeof(int),                     &M,         VALUE,
        sizeof(int),                     &N,         VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), INOUT,
        0);
}
