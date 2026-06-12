/**
 *
 * @file quark/codelet_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm Quark codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_ztrmm_quark(Quark *quark)
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    quark_unpack_args_9(quark, side, uplo, transA, diag, M, N, alpha, tileA, tileB);
    TCORE_ztrmm(side, uplo,
        transA, diag,
        M, N,
        alpha, tileA,
        tileB);
}

void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlaset( options, ChamUpperLower, m, n,
                            alpha, alpha, B, Bm, Bn );
        return;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TRMM;
    QUARK_Insert_Task(opt->quark, CORE_ztrmm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &side,      VALUE,
        sizeof(int),                &uplo,      VALUE,
        sizeof(int),                &transA,    VALUE,
        sizeof(int),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),                 INOUT,
        0);

    (void)nb;
}
