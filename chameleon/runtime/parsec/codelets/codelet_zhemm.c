/**
 *
 * @file parsec/codelet_zhemm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm PaRSEC codelet
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
CORE_zhemm_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int LDC;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &M, &N, &alpha, &A, &LDA, &B, &LDB, &beta, &C, &LDC );

    CORE_zhemm( side, uplo, M, N,
                alpha, A, LDA,
                       B, LDB,
                beta,  C, LDC );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void
INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                   cham_side_t side, cham_uplo_t uplo,
                   int m, int n, int nb,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *B, int Bm, int Bn,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zhemm_parsec, options->priority, "hemm",
        sizeof(cham_side_t),                &side,    VALUE,
        sizeof(cham_uplo_t),                &uplo,    VALUE,
        sizeof(int),                       &m,       VALUE,
        sizeof(int),                       &n,       VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,   VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INPUT,
        sizeof(int), &(tileB->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t),         &beta,    VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}

void
INSERT_TASK_zhemm_Astat( const RUNTIME_option_t *options,
                         cham_side_t side, cham_uplo_t uplo,
                         int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    INSERT_TASK_zhemm( options, side, uplo, m, n, nb,
                       alpha, A, Am, An, B, Bm, Bn,
                       beta, C, Cm, Cn );
}
