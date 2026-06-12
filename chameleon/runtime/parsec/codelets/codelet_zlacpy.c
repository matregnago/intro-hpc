/**
 *
 * @file parsec/codelet_zlacpy.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlacpy_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int M;
    int N;
    const CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &A, &LDA, &B, &LDB );

    CORE_zlacpy( uplo, M, N, A, LDA, B, LDB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlacpy_parsec, options->priority, "lacpy",
        sizeof(cham_uplo_t), &uplo,        VALUE,
        sizeof(int),         &m,           VALUE,
        sizeof(int),         &n,           VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),         &(tileA->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | OUTPUT | AFFINITY,
        sizeof(int),         &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
}

static inline int
CORE_zlacpyx_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int M;
    int N;
    int displA;
    const CHAMELEON_Complex64_t *A;
    int LDA;
    int displB;
    CHAMELEON_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &displA, &A, &LDA, &displB, &B, &LDB );

    CORE_zlacpy( uplo, M, N, A + displA, LDA, B + displB, LDB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlacpyx_parsec, options->priority, "lacpyx",
        sizeof(cham_uplo_t), &uplo,   VALUE,
        sizeof(int),         &m,      VALUE,
        sizeof(int),         &n,      VALUE,
        sizeof(int),         &displA, VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),         &lda,    VALUE,
        sizeof(int),         &displB, VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | OUTPUT | AFFINITY,
        sizeof(int),         &ldb,    VALUE,
        PARSEC_DTD_ARG_END );
}
