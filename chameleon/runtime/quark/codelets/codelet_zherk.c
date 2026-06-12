/**
 *
 * @file quark/codelet_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zherk_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    double alpha;
    CHAM_tile_t *tileA;
    double beta;
    CHAM_tile_t *tileC;

    quark_unpack_args_8(quark, uplo, trans, n, k, alpha, tileA, beta, tileC);
    TCORE_zherk(uplo, trans,
               n, k,
               alpha, tileA,
               beta, tileC);
}

void INSERT_TASK_zherk(const RUNTIME_option_t *options,
                      cham_uplo_t uplo, cham_trans_t trans,
                      int n, int k, int nb,
                      double alpha, const CHAM_desc_t *A, int Am, int An,
                      double beta, const CHAM_desc_t *C, int Cm, int Cn)
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessC = ( beta == 0. ) ? OUTPUT : INOUT;

    DAG_CORE_HERK;
    QUARK_Insert_Task(opt->quark, CORE_zherk_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,      VALUE,
        sizeof(int),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),                     &alpha,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(double),                     &beta,      VALUE,
        sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),                 accessC,
        0);
}
