/**
 *
 * @file parsec/codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_blas PaRSEC wrapper
 *
 * @version 1.4.0
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlatro_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int M;
    int N;
    const CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &uplo, &trans, &M, &N, &A, &LDA, &B, &LDB);

    CORE_zlatro( uplo, trans, M, N,
                A, LDA, B, LDB);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlatro(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int m, int n, int mb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlatro_parsec, options->priority, "latro",
        sizeof(cham_uplo_t), &uplo,  VALUE,
        sizeof(cham_trans_t), &trans, VALUE,
        sizeof(int),        &m,     VALUE,
        sizeof(int),        &n,     VALUE,
        PASSED_BY_REF,       RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), OUTPUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );

    (void)mb;
}
