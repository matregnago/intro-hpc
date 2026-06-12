/**
 *
 * @file quark/codelet_zaxpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy Quark codelet
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

void CORE_zaxpy_quark(Quark *quark)
{
    int M;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    int incA;
    CHAM_tile_t *tileB;
    int incB;

    quark_unpack_args_6( quark, M, alpha, tileA, incA, tileB, incB );
    TCORE_zaxpy(M, alpha, tileA, incA, tileB, incB);
}

void INSERT_TASK_zaxpy(const RUNTIME_option_t *options,
                      int M, CHAMELEON_Complex64_t alpha,
                      const CHAM_desc_t *A, int Am, int An, int incA,
                      const CHAM_desc_t *B, int Bm, int Bn, int incB)
{
    if ( alpha == 0. ) {
        return;
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_AXPY;
    QUARK_Insert_Task(opt->quark, CORE_zaxpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &M,         VALUE,
        sizeof(CHAMELEON_Complex64_t),      &alpha,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                        &incA,      VALUE,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), INOUT,
        sizeof(int),                        &incB,      VALUE,
        0);
}
