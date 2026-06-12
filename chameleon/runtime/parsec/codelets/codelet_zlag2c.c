/**
 *
 * @file parsec/codelet_zlag2c.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_zc.h"
#include "coreblas/coreblas_zc.h"

static inline int
CORE_zlag2c_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int info;
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex32_t *B;
    int ldb;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &B, &ldb );

    CORE_zlag2c( m, n, A, lda, B, ldb, &info );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlag2c_parsec, options->priority, "lag2c",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, CHAMELEON_Complex32_t, Bm, Bn ),     OUTPUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
}

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
static inline int
CORE_clag2z_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    CHAMELEON_Complex32_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &B, &ldb );

    CORE_clag2z( m, n, A, lda, B, ldb );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_clag2z_parsec, options->priority, "lag2z",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex32_t, Am, An ),     INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | OUTPUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
}
