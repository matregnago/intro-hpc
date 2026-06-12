/**
 *
 * @file parsec/codelet_ztrtri.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri PaRSEC codelet
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
CORE_ztrtri_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    cham_diag_t diag;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    parsec_dtd_unpack_args(
        this_task, &uplo, &diag, &N, &A, &LDA, &iinfo, &sequence, &request );

    CORE_ztrtri( uplo, diag, N, A, LDA, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_diag_t diag,
                        int n, int nb,
                        const CHAM_desc_t *A, int Am, int An,
                        int iinfo )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztrtri_parsec, options->priority, "trtri",
        sizeof(cham_uplo_t),                 &uplo,                  VALUE,
        sizeof(cham_diag_t),                 &diag,                  VALUE,
        sizeof(int),                 &n,                     VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),                 &iinfo,                 VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),   VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),    VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
