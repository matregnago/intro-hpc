/**
 *
 * @file parsec/codelet_zsyssq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyssq PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zsyssq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    parsec_dtd_unpack_args(
        this_task, &storev, &uplo, &n, &A, &lda, &SCALESUMSQ );

    CORE_zsyssq( storev, uplo, n, A, lda, SCALESUMSQ );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                        cham_store_t storev, cham_uplo_t uplo, int n,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zsyssq_parsec, options->priority, "syssq",
        sizeof(cham_store_t),   &storev,                VALUE,
        sizeof(cham_uplo_t),            &uplo,                  VALUE,
        sizeof(int),            &n,                     VALUE,
        PASSED_BY_REF,          RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,          RTBLKADDR( SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn ), chameleon_parsec_get_arena_index( SCALESUMSQ ) | INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
