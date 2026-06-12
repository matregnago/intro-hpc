/**
 *
 * @file quark/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy Quark codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static inline void CORE_zlacpy_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M, N;
    CHAM_tile_t *tileA, *tileB;

    quark_unpack_args_5(quark, uplo, M, N, tileA, tileB);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    TCORE_zlacpy( uplo, M, N, tileA, tileB );
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zlacpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),   &uplo,   VALUE,
        sizeof(int),   &m,      VALUE,
        sizeof(int),   &n,      VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             OUTPUT,
        0);
}

static inline void CORE_zlacpyx_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M, N;
    int displA, displB;
    int LDA, LDB;
    CHAM_tile_t *tileA, *tileB;

    quark_unpack_args_9(quark, uplo, M, N, displA, tileA, LDA, displB, tileB, LDB);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    TCORE_zlacpyx( uplo, M, N, displA, tileA, LDA, displB, tileB, LDB );
}

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zlacpyx_quark, (Quark_Task_Flags*)opt,
        sizeof(int),   &uplo,   VALUE,
        sizeof(int),   &m,      VALUE,
        sizeof(int),   &n,      VALUE,
        sizeof(int),   &displA, VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(int),   &lda,    VALUE,
        sizeof(int),   &displB, VALUE,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             OUTPUT,
        sizeof(int),   &ldb,    VALUE,
        0);
}
