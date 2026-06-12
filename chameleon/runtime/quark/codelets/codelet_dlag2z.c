/**
 *
 * @file quark/codelet_dlag2z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlag2z Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static inline void
CORE_dlag2z_quark( Quark *quark )
{
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    quark_unpack_args_5( quark, uplo, m, n, tileA, tileB );
    TCORE_dlag2z( uplo, m, n, tileA, tileB );
}

void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    //DAG_CORE_DLAG2Z;
    QUARK_Insert_Task(
        opt->quark, CORE_dlag2z_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_uplo_t), &uplo,    VALUE,
        sizeof(int),         &m,         VALUE,
        sizeof(int),         &n,         VALUE,
        sizeof(void*), RTBLKADDR(A, double, Am, An),                INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), OUTPUT,
        0);
}
