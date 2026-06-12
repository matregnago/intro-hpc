/**
 *
 * @file quark/codelet_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange Quark codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zlange_quark(Quark *quark)
{
    CHAM_tile_t *tileNorm;
    cham_normtype_t norm;
    int M;
    int N;
    CHAM_tile_t *tileA;
    double *work;

    quark_unpack_args_6( quark, norm, M, N, tileA, work, tileNorm );
    TCORE_zlange( norm, M, N, tileA, work, tileNorm->mat );
}

void INSERT_TASK_zlange(const RUNTIME_option_t *options,
                       cham_normtype_t norm, int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANGE;
    int szeW = chameleon_max( m, n );
    QUARK_Insert_Task(
        opt->quark, CORE_zlange_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &norm,  VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);

    (void)nb;
}

void CORE_zlange_max_quark(Quark *quark)
{
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileNorm;
    double *A, *norm;

    quark_unpack_args_2(quark, tileA, tileNorm);
    A = tileA->mat;
    norm = tileNorm->mat;

    if ( A[0] > *norm ) {
        *norm = A[0];
    }
}

void INSERT_TASK_zlange_max(const RUNTIME_option_t *options,
                           const CHAM_desc_t *A, int Am, int An,
                           const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANGE_MAX;
    QUARK_Insert_Task(
        opt->quark, CORE_zlange_max_quark, (Quark_Task_Flags*)opt,
        sizeof(void*), RTBLKADDR(A, double, Am, An), INPUT,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);
}

