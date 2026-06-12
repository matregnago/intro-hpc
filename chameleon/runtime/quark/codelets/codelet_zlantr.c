/**
 *
 * @file quark/codelet_zlantr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlantr Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zlantr_quark(Quark *quark)
{
    CHAM_tile_t *tileNorm;
    cham_normtype_t norm;
    cham_uplo_t     uplo;
    cham_diag_t     diag;
    int M;
    int N;
    CHAM_tile_t *tileA;
    double *work;

    quark_unpack_args_8(quark, norm, uplo, diag, M, N, tileA, work, tileNorm );
    TCORE_zlantr( norm, uplo, diag, M, N, tileA, work, tileNorm->mat );
}

void INSERT_TASK_zlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANTR;
    int szeW = chameleon_max( 1, n );
    QUARK_Insert_Task(
        opt->quark, CORE_zlantr_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &norm,  VALUE,
        sizeof(int),              &uplo,  VALUE,
        sizeof(int),              &diag,  VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);

    (void)nb;
}
