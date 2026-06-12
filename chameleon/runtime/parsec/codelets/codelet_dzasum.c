/**
 *
 * @file parsec/codelet_dzasum.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dzasum PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_dzasum_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *work;

    parsec_dtd_unpack_args(
        this_task, &storev, &uplo, &M, &N, &A, &lda, &work );

    CORE_dzasum( storev, uplo, M, N, A, lda, work );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_dzasum(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, int M, int N,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_dzasum_parsec, options->priority, "dzasum",
        sizeof(int),    &storev,                           VALUE,
        sizeof(int),    &uplo,                             VALUE,
        sizeof(int),           &M,                                VALUE,
        sizeof(int),           &N,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),     INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
