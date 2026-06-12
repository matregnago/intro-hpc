/**
 *
 * @file parsec/codelet_dlag2z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlag2z PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_dlag2z_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int m;
    int n;
    double *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;

    parsec_dtd_unpack_args( this_task, &uplo, &m, &n, &A, &lda, &B, &ldb );
    CORE_dlag2z( uplo, m, n, A, lda, B, ldb );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_dlag2z_parsec, options->priority, "dlag2z",
        sizeof(cham_uplo_t), &uplo,       VALUE,
        sizeof(int),         &m,          VALUE,
        sizeof(int),         &n,          VALUE,
        PASSED_BY_REF,        RTBLKADDR( A, double, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),        &(tileA->ld), VALUE,
        PASSED_BY_REF,        RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | OUTPUT,
        sizeof(int),        &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
}
