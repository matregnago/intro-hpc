/**
 *
 * @file parsec/codelet_zgeadd.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgeadd_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_trans_t trans;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &trans, &M, &N, &alpha, &A, &LDA, &beta, &B, &LDB );

    CORE_zgeadd( trans, M, N, alpha, A, LDA, beta, B, LDB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

/**
 ******************************************************************************
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 * @brief Adds two general matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = ChamNoTrans:   op(A) = A
 *          = ChamTrans:     op(A) = A'
 *          = ChamConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = ChamNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = ChamNoTrans, and k=N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgeadd_parsec, options->priority, "geadd",
        sizeof(cham_trans_t),        &trans, VALUE,
        sizeof(int),               &m,     VALUE,
        sizeof(int),               &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,              RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t), &beta,  VALUE,
        PASSED_BY_REF,              RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
