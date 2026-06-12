/**
 *
 * @file core_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zherfb CPU kernel
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/lapacke.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zherfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q^T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zgeqrt. Only ChamLower supported!
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         TODO comment
 *
 * @param[in] n
 *         The number of rows/columns of the tile C.  N >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q. K >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in] A
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgeqrt in the first k columns of its array argument A.
 *
 * @param[in] lda
 *         The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the symmetric N-by-N tile C.
 *         On exit, C is overwritten by Q^T*C*Q.
 *
 * @param[in] ldc
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] ldwork
 *         The dimension of the array WORK. LDWORK >= max(1,N);
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
#if defined(CHAMELEON_HAVE_WEAK)
#pragma weak CORE_zherfb = PCORE_zherfb
#define CORE_zherfb PCORE_zherfb
#define CORE_zunmlq PCORE_zunmlq
#define CORE_zunmqr PCORE_zunmqr
int  CORE_zunmlq(cham_side_t side, cham_trans_t trans,
                 int M, int N, int IB, int K,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *C, int LDC,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(cham_side_t side, cham_trans_t trans,
                 int M, int N, int K, int IB,
                 const CHAMELEON_Complex64_t *V, int LDV,
                 const CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *C, int LDC,
                 CHAMELEON_Complex64_t *WORK, int LDWORK);
#endif
int CORE_zherfb( cham_uplo_t uplo, int n,
                 int k, int ib, int nb,
                 const CHAMELEON_Complex64_t *A, int lda,
                 const CHAMELEON_Complex64_t *T, int ldt,
                 CHAMELEON_Complex64_t *C, int ldc,
                 CHAMELEON_Complex64_t *WORK, int ldwork )
{
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        coreblas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        coreblas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        coreblas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        coreblas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < chameleon_max(1,n)) && (n > 0) ) {
        coreblas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < chameleon_max(1,ib)) && (ib > 0) ) {
        coreblas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < chameleon_max(1,n)) && (n > 0) ) {
        coreblas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == ChamLower) {
        /* Left */
        CORE_zunmqr(ChamLeft, ChamConjTrans, n, n, k, ib,
                    A, lda, T, ldt, C, ldc, WORK, ldwork);
        /* Right */
        CORE_zunmqr(ChamRight, ChamNoTrans, n, n, k, ib,
                    A, lda, T, ldt, C, ldc, WORK, ldwork);
    }
    else {
        /* Right */
        CORE_zunmlq(ChamRight, ChamConjTrans, n, n, k, ib,
                    A, lda, T, ldt, C, ldc, WORK, ldwork);
        /* Left */
        CORE_zunmlq(ChamLeft, ChamNoTrans, n, n, k, ib,
                    A, lda, T, ldt, C, ldc, WORK, ldwork);
    }
    return 0;
}
