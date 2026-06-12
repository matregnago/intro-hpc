/**
 *
 * @file quark/codelet_zgesum.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesum Quark codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgesum_quark(Quark *quark)
{
    cham_store_t storev;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    quark_unpack_args_5( quark, storev, m, n, tileA, tileW );
    TCORE_zgesum( storev, m, n, tileA, tileW );
}

void INSERT_TASK_zgesum( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SUMS, int SUMSm, int SUMSn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GESUM;
    QUARK_Insert_Task(opt->quark, CORE_zgesum_quark, (Quark_Task_Flags*)opt,
                      sizeof(cham_store_t),            &storev, VALUE,
                      sizeof(int),                     &m,      VALUE,
                      sizeof(int),                     &n,      VALUE,
                      sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
                      sizeof(void*), RTBLKADDR(SUMS, CHAMELEON_Complex64_t, SUMSm, SUMSn), INOUT,
                      0);
}
