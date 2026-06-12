/**
 *
 * @file parsec/codelet_zplssq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq PaRSEC codelet
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
CORE_zplssq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_store_t storev;
    int M;
    int N;
    double *SCLSSQ_IN;
    double *SCLSSQ_OUT;

    parsec_dtd_unpack_args(
        this_task, &storev, &M, &N, &SCLSSQ_IN, &SCLSSQ_OUT );

    CORE_zplssq(storev, M, N, SCLSSQ_IN, SCLSSQ_OUT);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplssq_parsec, options->priority, "plssq",
        sizeof(int),           &storev,                           VALUE,
        sizeof(int),           &M,                                VALUE,
        sizeof(int),           &N,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn ), chameleon_parsec_get_arena_index( SCALESUMSQ) | INPUT,
        PASSED_BY_REF,         RTBLKADDR( SCLSSQ, double, SCLSSQm, SCLSSQn ), chameleon_parsec_get_arena_index( SCLSSQ) | INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}

static inline int
CORE_zplssq2_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    int N;
    double *RESULT;

    parsec_dtd_unpack_args(
        this_task, &N, &RESULT );

    CORE_zplssq2(N, RESULT);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplssq2_parsec, options->priority, "plssq2",
        sizeof(int),           &N,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( RESULT, double, RESULTm, RESULTn ), chameleon_parsec_get_arena_index( RESULT) | INOUT | AFFINITY,
        PARSEC_DTD_ARG_END );
}
