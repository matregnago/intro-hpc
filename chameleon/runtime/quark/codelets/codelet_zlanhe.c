/**
 *
 * @file quark/codelet_zlanhe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlanhe Quark codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zlanhe_quark(Quark *quark)
{
    CHAM_tile_t *tileNorm;
    cham_normtype_t norm;
    cham_uplo_t uplo;
    int N;
    CHAM_tile_t *tileA;
    double *work;

    quark_unpack_args_6(quark, norm, uplo, N, tileA, work, tileNorm );
    TCORE_zlanhe( norm, uplo, N, tileA, work, tileNorm->mat );
}

void INSERT_TASK_zlanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANHE;
    int szeW = chameleon_max( 1, n );
    QUARK_Insert_Task(
        opt->quark, CORE_zlanhe_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &norm,  VALUE,
        sizeof(int),              &uplo,  VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(void*), RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);

    (void)nb;
}
