/**
 *
 * @file quark/codelet_zgemv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemv Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgemv_quark(Quark *quark)
{
    cham_trans_t trans;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileX;
    int incX;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileY;
    int incY;

    quark_unpack_args_10( quark, trans, m, n, alpha, tileA, tileX, incX, beta, tileY, incY );
    TCORE_zgemv( trans, m, n,
                 alpha, tileA, tileX, incX,
                 beta,  tileY, incY );
}

void INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessY = ( beta == 0. ) ? OUTPUT : INOUT;

    /* DAG_CORE_GEMV; */
    QUARK_Insert_Task(
        opt->quark, CORE_zgemv_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_trans_t),          &trans, VALUE,
        sizeof(int),                   &m,     VALUE,
        sizeof(int),                   &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t), &alpha, VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(void*), RTBLKADDR(X, CHAMELEON_Complex64_t, Xm, Xn), INPUT,
        sizeof(int),                   &incX,  VALUE,
        sizeof(CHAMELEON_Complex64_t), &beta,  VALUE,
        sizeof(void*), RTBLKADDR(Y, CHAMELEON_Complex64_t, Ym, Yn), accessY,
        sizeof(int),                   &incY,  VALUE,
        0);
}
