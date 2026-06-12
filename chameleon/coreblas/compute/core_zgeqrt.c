/**
 *
 * @file core_zgeqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgeqrt CPU kernel
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
 *  CORE_zgeqrt computes a QR factorization of a complex M-by-N tile A:
 *  A = Q * R.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the min(M,N)-by-N upper trapezoidal tile R (R is
 *         upper triangular if M >= N); the elements below the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zgeqrt(int M, int N, int IB,
                CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *T, int LDT,
                CHAMELEON_Complex64_t *TAU,
                CHAMELEON_Complex64_t *WORK)
{
    int i, k, sb;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ((LDT < chameleon_max(1,IB)) && (IB > 0)) {
        coreblas_error(7, "Illegal value of LDT");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    k = chameleon_min(M, N);

    for(i = 0; i < k; i += IB) {
        sb = chameleon_min(IB, k-i);

        LAPACKE_zgeqr2_work(LAPACK_COL_MAJOR, M-i, sb,
                            &A[LDA*i+i], LDA, &TAU[i], WORK);

        LAPACKE_zlarft_work(LAPACK_COL_MAJOR,
            chameleon_lapack_const(ChamDirForward),
            chameleon_lapack_const(ChamColumnwise),
            M-i, sb,
            &A[LDA*i+i], LDA, &TAU[i],
            &T[LDT*i], LDT);

        if (N > i+sb) {
            LAPACKE_zlarfb_work(
                LAPACK_COL_MAJOR,
                chameleon_lapack_const(ChamLeft),
                chameleon_lapack_const(ChamConjTrans),
                chameleon_lapack_const(ChamDirForward),
                chameleon_lapack_const(ChamColumnwise),
                M-i, N-i-sb, sb,
                &A[LDA*i+i],      LDA,
                &T[LDT*i],        LDT,
                &A[LDA*(i+sb)+i], LDA,
                WORK, N-i-sb);
        }
    }
    return CHAMELEON_SUCCESS;
}


