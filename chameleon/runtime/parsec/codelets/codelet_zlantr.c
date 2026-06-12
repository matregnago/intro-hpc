/**
 *
 * @file parsec/codelet_zlantr.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlantr PaRSEC codelet
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
CORE_zlantr_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_normtype_t norm;
    cham_uplo_t uplo;
    cham_diag_t diag;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;
    double *work;
    double *normA;

    parsec_dtd_unpack_args(
        this_task, &norm, &uplo, &diag, &M, &N, &A, &LDA, &work, &normA );

    CORE_zlantr( norm, uplo, diag, M, N, A, LDA, work, normA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlantr(const RUNTIME_option_t *options,
                       cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                       int M, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    int szeW = chameleon_max( 1, N );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlantr_parsec, options->priority, "lantr",
        sizeof(cham_normtype_t), &norm,       VALUE,
        sizeof(cham_uplo_t), &uplo,       VALUE,
        sizeof(cham_diag_t), &diag,       VALUE,
        sizeof(int), &M,          VALUE,
        sizeof(int), &N,          VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(double)*szeW,           NULL,           SCRATCH,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),            OUTPUT | AFFINITY,
        PARSEC_DTD_ARG_END );

    (void)NB;
}
