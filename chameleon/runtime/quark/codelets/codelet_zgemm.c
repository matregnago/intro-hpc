/**
 *
 * @file quark/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
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

void CORE_zgemm_quark(Quark *quark)
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileC;

    quark_unpack_args_10(quark, transA, transB, m, n, k, alpha, tileA, tileB, beta, tileC);
    TCORE_zgemm( transA, transB,
                 m, n, k,
                 alpha, tileA, tileB,
                 beta,  tileC );
}

void
INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                   cham_trans_t transA, cham_trans_t transB,
                   int m, int n, int k, int nb,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *B, int Bm, int Bn,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessC = ( beta == 0. ) ? OUTPUT : INOUT;

    DAG_CORE_GEMM;
    QUARK_Insert_Task(opt->quark, CORE_zgemm_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                &transA,    VALUE,
                      sizeof(int),                &transB,    VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(int),                        &k,         VALUE,
                      sizeof(CHAMELEON_Complex64_t),         &alpha,     VALUE,
                      sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
                      sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),                 INPUT,
                      sizeof(CHAMELEON_Complex64_t),         &beta,      VALUE,
                      sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),                 accessC,
                      0);
}

void
INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    INSERT_TASK_zgemm( options, transA, transB, m, n, k, nb,
                       alpha, A, Am, An, B, Bm, Bn,
                       beta, C, Cm, Cn );
}
