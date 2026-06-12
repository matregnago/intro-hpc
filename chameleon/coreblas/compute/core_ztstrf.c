/**
 *
 * @file core_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztstrf CPU kernel
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
 * @date 2025-01-29
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include <math.h>

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_ztstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular NB-by-N tile U on top of a M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] NB
 *
 * @param[in,out] U
 *         On entry, the NB-by-N upper triangular tile.
 *         On exit, the new factor U from the factorization
 *
 * @param[in] LDU
 *         The leading dimension of the array U.  LDU >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factor L from the factorization
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[in,out] L
 *         On entry, the IB-by-N lower triangular tile.
 *         On exit, the interchanged rows form the tile A in case of pivoting.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,IB).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile U was interchanged with row IPIV(i) of the tile A.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *
 * @param[out] INFO
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 * @retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 */


int CORE_ztstrf(int M, int N, int IB, int NB,
                CHAMELEON_Complex64_t *U, int LDU,
                CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *L, int LDL,
                int *IPIV,
                CHAMELEON_Complex64_t *WORK, int LDWORK,
                int *INFO)
{
    static CHAMELEON_Complex64_t zzero = 0.0;
    static CHAMELEON_Complex64_t mzone =-1.0;

    CHAMELEON_Complex64_t alpha;
    int i, j, ii, sb;
    int im, ip;

    /* Check input arguments */
    *INFO = 0;
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
    if ((LDU < chameleon_max(1,NB)) && (NB > 0)) {
        coreblas_error(6, "Illegal value of LDU");
        return -6;
    }
    if ((LDA < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if ((LDL < chameleon_max(1,IB)) && (IB > 0)) {
        coreblas_error(10, "Illegal value of LDL");
        return -10;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    /* Set L to 0 */
    memset( L, 0, sizeof(CHAMELEON_Complex64_t) * LDL * N );

    ip = 0;
    for (ii = 0; ii < N; ii += IB) {
        sb = chameleon_min(N-ii, IB);

        for (i = 0; i < sb; i++) {
            im = cblas_izamax(M, &A[LDA*(ii+i)], 1);
            IPIV[ip] = ii+i+1;

            if (cabs(A[LDA*(ii+i)+im]) > cabs(U[LDU*(ii+i)+ii+i])) {
                /*
                 * Swap behind.
                 */
                cblas_zswap(i, &L[LDL*ii+i], LDL, &WORK[im], LDWORK );
                /*
                 * Swap ahead.
                 */
                cblas_zswap(sb-i, &U[LDU*(ii+i)+ii+i], LDU, &A[LDA*(ii+i)+im], LDA );
                /*
                 * Set IPIV.
                 */
                IPIV[ip] = NB + im + 1;

                for (j = 0; j < i; j++) {
                    A[LDA*(ii+j)+im] = zzero;
                }
            }

            if ((*INFO == 0) && (cabs(A[LDA*(ii+i)+im]) == zzero)
                && (cabs(U[LDU*(ii+i)+ii+i]) == zzero)) {
                *INFO = ii+i+1;
            }

            alpha = ((CHAMELEON_Complex64_t)1. / U[LDU*(ii+i)+ii+i]);
            cblas_zscal(M, CBLAS_SADDR(alpha), &A[LDA*(ii+i)], 1);
            cblas_zcopy(M, &A[LDA*(ii+i)], 1, &WORK[LDWORK*i], 1);
            cblas_zgeru(
                CblasColMajor, M, sb-i-1,
                CBLAS_SADDR(mzone), &A[LDA*(ii+i)], 1,
                &U[LDU*(ii+i+1)+ii+i], LDU,
                &A[LDA*(ii+i+1)], LDA );
            ip = ip+1;
        }
        /*
         * Apply the subpanel to the rest of the panel.
         */
        if(ii+i < N) {
            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] - ii;
                }
            }

            CORE_zssssm(
                NB, N-(ii+sb), M, N-(ii+sb), sb, sb,
                &U[LDU*(ii+sb)+ii], LDU,
                &A[LDA*(ii+sb)], LDA,
                &L[LDL*ii], LDL,
                WORK, LDWORK, &IPIV[ii]);

            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] + ii;
                }
            }
        }
    }
    return CHAMELEON_SUCCESS;
}


