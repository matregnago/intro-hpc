/**
 *
 * @file parsec/codelet_zgemv.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemv PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgemv_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    cham_trans_t trans;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *X;
    int incX;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *Y;
    int incY;

    parsec_dtd_unpack_args(
        this_task, &trans, &m, &n, &alpha, &A, &lda, &X, &incX, &beta, &Y, &incY );

    CORE_zgemv( trans, m, n,
                alpha, A, lda,
                       X, incX,
                beta,  Y, incY );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void
INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                   cham_trans_t trans, int m, int n,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *X, int Xm, int Xn, int incX,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgemv_parsec, options->priority, "zgemv",
        sizeof(cham_trans_t),          &trans,       VALUE,
        sizeof(int),                   &m,           VALUE,
        sizeof(int),                   &n,           VALUE,
        sizeof(CHAMELEON_Complex64_t), &alpha,       VALUE,
        PASSED_BY_REF, RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),                   &(tileA->ld), VALUE,
        PASSED_BY_REF, RTBLKADDR( X, CHAMELEON_Complex64_t, Xm, Xn ), chameleon_parsec_get_arena_index( X ) | INPUT,
        sizeof(int),                   &incX,        VALUE,
        sizeof(CHAMELEON_Complex64_t), &beta,        VALUE,
        PASSED_BY_REF, RTBLKADDR( Y, CHAMELEON_Complex64_t, Ym, Yn ), chameleon_parsec_get_arena_index( Y ) | INOUT | AFFINITY,
        sizeof(int),                   &incY,        VALUE,
        PARSEC_DTD_ARG_END );
}
