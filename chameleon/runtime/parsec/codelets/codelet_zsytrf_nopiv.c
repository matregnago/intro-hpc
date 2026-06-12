/**
 *
 * @file parsec/codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zsytrf_nopiv_parsec( parsec_execution_stream_t *context,
                          parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    int iinfo;

    parsec_dtd_unpack_args(
        this_task, &uplo, &n, &A, &lda, &iinfo );

    CORE_zsytf2_nopiv( uplo, n, A, lda );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zsytrf_nopiv(const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zsytrf_nopiv_parsec, options->priority, "sytrf_nopiv",
        sizeof(cham_uplo_t),              &uplo,                VALUE,
        sizeof(int),                     &n,                   VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),                     &iinfo,               VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
