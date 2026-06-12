/**
 *
 * @file core_zttqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zttqrt CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Dulceneia Becker
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
 *  CORE_zttqrt computes a QR factorization of a rectangular matrix
 *  formed by coupling a complex N-by-N upper triangular tile A1
 *  on top of a complex M-by-N upper trapezoidal tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
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
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A2(1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A2.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A1 and A2.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N upper triangular tile A2.
 *         On exit, the elements on and above the diagonal of the array
 *         with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M).
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
 * @param[in,out] WORK
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zttqrt(int M, int N, int IB,
                CHAMELEON_Complex64_t *A1, int LDA1,
                CHAMELEON_Complex64_t *A2, int LDA2,
                CHAMELEON_Complex64_t *T, int LDT,
                CHAMELEON_Complex64_t *TAU, CHAMELEON_Complex64_t *WORK)
{
    static int                ione  = 1;
    static CHAMELEON_Complex64_t zone  = 1.0;
    static CHAMELEON_Complex64_t zzero = 0.0;

    CHAMELEON_Complex64_t alpha;
    int i, j, l, ii, sb, mi, ni;

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
        coreblas_error(7, "Illegal value of LDA2");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    /* TODO: Need to check why some cases require
     *  this to not have uninitialized values */
    CORE_zlaset( ChamUpperLower, IB, N,
                 0., 0., T, LDT);

    for (ii = 0; ii < N; ii += IB) {
        sb = chameleon_min(N-ii, IB);
        for (i = 0; i < sb; i++) {
            j  = ii + i;
            mi = chameleon_min( j + 1, M );
            ni = sb-i-1;

            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:mi, II*IB+I ).
             */
            LAPACKE_zlarfg_work(
                    mi+1, &A1[LDA1*j+j], &A2[LDA2*j], ione, &TAU[j]);

            if (ni > 0) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left.
                 */
                cblas_zcopy(
                    ni,
                    &A1[LDA1*(j+1)+j], LDA1,
                    WORK, 1);

#if defined(PRECISION_z) || defined(PRECISION_c)
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                cblas_zgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)ChamConjTrans,
                    mi, ni,
                    CBLAS_SADDR(zone), &A2[LDA2*(j+1)], LDA2,
                                       &A2[LDA2*j],     1,
                    CBLAS_SADDR(zone), WORK,            1);
#if defined(PRECISION_z) || defined(PRECISION_c)
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                alpha = -conj(TAU[j]);
                cblas_zaxpy(
                    ni, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(j+1)+j], LDA1);
#if defined(PRECISION_z) || defined(PRECISION_c)
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                cblas_zgerc(
                    CblasColMajor, mi, ni,
                    CBLAS_SADDR(alpha), &A2[LDA2*j], 1,
                    WORK, 1,
                    &A2[LDA2*(j+1)], LDA2);
            }

            /*
             * Calculate T
             *
             * T(0:i-1, j) = alpha * A2(0:M-1, ii:j-1)' * A2(0:M-1, j)
             */

            if ( i > 0 ) {

                l = chameleon_min(i, chameleon_max(0, M-ii));
                alpha = -(TAU[j]);

                CORE_zpemv(
                        ChamConjTrans, ChamColumnwise,
                        chameleon_min(j, M), i, l,
                        alpha, &A2[LDA2*ii], LDA2,
                               &A2[LDA2*j],  1,
                        zzero, &T[LDT*j],    1,
                        WORK);

                /* T(0:i-1, j) = T(0:i-1, ii:j-1) * T(0:i-1, j) */
                cblas_ztrmv(
                        CblasColMajor, (CBLAS_UPLO)ChamUpper,
                        (CBLAS_TRANSPOSE)ChamNoTrans,
                        (CBLAS_DIAG)ChamNonUnit,
                        i, &T[LDT*ii], LDT,
                           &T[LDT*j], 1);

            }

            T[LDT*j+i] = TAU[j];
        }

        /* Apply Q' to the rest of the matrix to the left  */
        if (N > ii+sb) {
            mi = chameleon_min(ii+sb, M);
            ni = N-(ii+sb);
            l  = chameleon_min(sb, chameleon_max(0, mi-ii));
            CORE_zparfb(
                ChamLeft, ChamConjTrans,
                ChamDirForward, ChamColumnwise,
                IB, ni, mi, ni, sb, l,             //replaced sb by IB
                &A1[LDA1*(ii+sb)+ii], LDA1,
                &A2[LDA2*(ii+sb)], LDA2,
                &A2[LDA2*ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, sb);
        }
    }
    return CHAMELEON_SUCCESS;
}


