/**
 *
 * @file core_ztslqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztslqt CPU kernel
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
 *  CORE_ztslqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M
 *  lower triangular tile A1 and a complex M-by-N tile A2:
 *
 *    | A1 A2 | = L * Q
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, all the elements with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
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

int CORE_ztslqt(int M, int N, int IB,
                CHAMELEON_Complex64_t *A1, int LDA1,
                CHAMELEON_Complex64_t *A2, int LDA2,
                CHAMELEON_Complex64_t *T, int LDT,
                CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK)
{
    static CHAMELEON_Complex64_t zone  = 1.0;
    static CHAMELEON_Complex64_t zzero = 0.0;

    CHAMELEON_Complex64_t alpha;
    int i, ii, sb;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA2 < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA2");
        return -8;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    for(ii = 0; ii < M; ii += IB) {
        sb = chameleon_min(M-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate A( II*IB+I, II*IB+I:N ).
             */
#if defined(PRECISION_z) || defined(PRECISION_c)
            LAPACKE_zlacgv_work(N, &A2[ii+i], LDA2);
            LAPACKE_zlacgv_work(1, &A1[LDA1*(ii+i)+ii+i], LDA1);
#endif
            LAPACKE_zlarfg_work(N+1, &A1[LDA1*(ii+i)+ii+i], &A2[ii+i], LDA2, &TAU[ii+i]);

            alpha = -(TAU[ii+i]);
            if (ii+i+1 < M) {
                /*
                 * Apply H( II+I-1 ) to A( II+I:II+IB-1, II+I-1:N  ) from the right.
                 */
                cblas_zcopy(
                    sb-i-1,
                    &A1[LDA1*(ii+i)+(ii+i+1)], 1,
                    WORK, 1);

                cblas_zgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)ChamNoTrans,
                    sb-i-1, N,
                    CBLAS_SADDR(zone), &A2[ii+i+1], LDA2,
                    &A2[ii+i], LDA2,
                    CBLAS_SADDR(zone), WORK, 1);

                cblas_zaxpy(
                    sb-i-1, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i)+ii+i+1], 1);

                cblas_zgerc(
                    CblasColMajor, sb-i-1, N,
                    CBLAS_SADDR(alpha), WORK, 1,
                    &A2[ii+i], LDA2,
                    &A2[ii+i+1], LDA2);
            }
            /*
             * Calculate T.
             */
            cblas_zgemv(
                CblasColMajor, (CBLAS_TRANSPOSE)ChamNoTrans, i, N,
                CBLAS_SADDR(alpha), &A2[ii], LDA2,
                &A2[ii+i], LDA2,
                CBLAS_SADDR(zzero), &T[LDT*(ii+i)], 1);
#if defined(PRECISION_z) || defined(PRECISION_c)
            LAPACKE_zlacgv_work(N, &A2[ii+i], LDA2 );
            LAPACKE_zlacgv_work(1, &A1[LDA1*(ii+i)+ii+i], LDA1 );
#endif
            cblas_ztrmv(
                CblasColMajor, (CBLAS_UPLO)ChamUpper,
                (CBLAS_TRANSPOSE)ChamNoTrans, (CBLAS_DIAG)ChamNonUnit, i,
                &T[LDT*ii], LDT,
                &T[LDT*(ii+i)], 1);

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }
        if (M > ii+sb) {
            CORE_ztsmlq(
                ChamRight, ChamConjTrans,
                M-(ii+sb), sb, M-(ii+sb), N, IB, IB,
                &A1[LDA1*ii+ii+sb], LDA1,
                &A2[ii+sb], LDA2,
                &A2[ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, M);
        }
    }
    return CHAMELEON_SUCCESS;
}


