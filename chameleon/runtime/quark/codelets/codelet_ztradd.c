/**
 *
 * @file quark/codelet_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_ztradd_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileB;

    quark_unpack_args_8(quark, uplo, trans, M, N, alpha, tileA, beta, tileB);
    TCORE_ztradd(uplo, trans, M, N, alpha, tileA, beta, tileB);
    return;
}

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessB = ( beta == 0. ) ? OUTPUT : INOUT;

    DAG_CORE_GEADD;
    QUARK_Insert_Task(opt->quark, CORE_ztradd_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                 &uplo,  VALUE,
        sizeof(int),                 &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha, VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),             INPUT,
        sizeof(CHAMELEON_Complex64_t),         &beta,   VALUE,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             accessB,
        0);

    (void)nb;
}

