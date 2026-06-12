/**
 *
 * @file parsec/codelet_zgesum.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesum PaRSEC codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgesum_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *SUMS;

    parsec_dtd_unpack_args(
        this_task, &storev, &m, &n, &A, &lda, &SUMS );

    CORE_zgesum( storev, m, n, A, lda, SUMS );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgesum( const RUNTIME_option_t *options,
                        cham_store_t storev, int m, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SUMS, int SUMSm, int SUMSn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgesum_parsec, options->priority, "gessum",
        sizeof(cham_store_t), &storev,       VALUE,
        sizeof(int),          &m,            VALUE,
        sizeof(int),          &n,            VALUE,
        PASSED_BY_REF,   RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,   RTBLKADDR( SUMS, CHAMELEON_Complex64_t, SUMSm, SUMSn ), chameleon_parsec_get_arena_index( SUMS ) | INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
