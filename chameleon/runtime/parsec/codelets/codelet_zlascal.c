/**
 *
 * @file parsec/codelet_zlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal PaRSEC codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlascal_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &alpha, &A, &LDA);

    CORE_zlascal( uplo, M, N, alpha, A, LDA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlascal(const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlascal_parsec, options->priority, "lascal",
        sizeof(cham_uplo_t),        &uplo,  VALUE,
        sizeof(int),               &m,     VALUE,
        sizeof(int),               &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,              RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}


