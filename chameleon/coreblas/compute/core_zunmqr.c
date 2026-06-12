/**
 *
 * @file core_zunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 */
 /**
 *
 * @brief Chameleon core_zunmqr CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zunmqr overwrites the general complex M-by-N tile C with
 *
 *                    SIDE = 'L'     SIDE = 'R'
 *    TRANS = 'N':      Q * C          C * Q
 *    TRANS = 'C':      Q^H * C       C * Q^H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zgeqrt. Q is of order M if SIDE = 'L' and of order N
 *  if SIDE = 'R'.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q^H from the Left;
 *         @arg ChamRight : apply Q or Q^H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  Transpose, apply Q^H.
 *
 * @param[in] M
 *         The number of rows of the tile C.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile C.  N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *         If SIDE = ChamLeft,  M >= K >= 0;
 *         if SIDE = ChamRight, N >= K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] A
 *         Dimension:  (LDA,K)
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgeqrt in the first k columns of its array argument A.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         If SIDE = ChamLeft,  LDA >= max(1,M);
 *         if SIDE = ChamRight, LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the M-by-N tile C.
 *         On exit, C is overwritten by Q*C or Q^T*C or C*Q^T or C*Q.
 *
 * @param[in] LDC
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *         If SIDE = ChamLeft,  LDWORK >= max(1,N);
 *         if SIDE = ChamRight, LDWORK >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zunmqr(cham_side_t side, cham_trans_t trans,
                int M, int N, int K, int IB,
                const CHAMELEON_Complex64_t *A, int LDA,
                const CHAMELEON_Complex64_t *T, int LDT,
                CHAMELEON_Complex64_t *C, int LDC,
                CHAMELEON_Complex64_t *WORK, int LDWORK)
{
    int i, kb;
    int i1, i3;
    int nq, nw;
    int ic = 0;
    int jc = 0;
    int ni = N;
    int mi = M;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        coreblas_error(1, "Illegal value of side");
        return -1;
    }
    /*
     * NQ is the order of Q and NW is the minimum dimension of WORK
     */
    if (side == ChamLeft) {
        nq = M;
        nw = N;
    }
    else {
        nq = N;
        nw = M;
    }

    if ((trans != ChamNoTrans) && (trans != ChamConjTrans)) {
        coreblas_error(2, "Illegal value of trans");
        return -2;
    }
    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ((K < 0) || (K > nq)) {
        coreblas_error(5, "Illegal value of K");
        return -5;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        coreblas_error(6, "Illegal value of IB");
        return -6;
    }
    if ((LDA < chameleon_max(1,nq)) && (nq > 0)) {
        coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if ((LDC < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(12, "Illegal value of LDC");
        return -12;
    }
    if ((LDWORK < chameleon_max(1,nw)) && (nw > 0)) {
        coreblas_error(14, "Illegal value of LDWORK");
        return -14;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0))
        return CHAMELEON_SUCCESS;

    if (((side == ChamLeft) && (trans != ChamNoTrans))
        || ((side == ChamRight) && (trans == ChamNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    for(i = i1; (i >- 1) && (i < K); i+=i3 ) {
        kb = chameleon_min(IB, K-i);

        if (side == ChamLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N - i;
            jc = i;
        }
        /*
         * Apply H or H'
         */
        LAPACKE_zlarfb_work(LAPACK_COL_MAJOR,
            chameleon_lapack_const(side),
            chameleon_lapack_const(trans),
            chameleon_lapack_const(ChamDirForward),
            chameleon_lapack_const(ChamColumnwise),
            mi, ni, kb,
            &A[LDA*i+i], LDA,
            &T[LDT*i], LDT,
            &C[LDC*jc+ic], LDC,
            WORK, LDWORK);
    }
    return CHAMELEON_SUCCESS;
}



