/**
 *
 * @file core_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgetrf_nopiv CPU kernel
 *
 * @version 1.4.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgetrf_nopiv computes an LU factorization of a general diagonal
 *  dominant M-by-N matrix A witout pivoting.
 *
 *  The factorization has the form
 *     A = L * U
 *  where L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 3 BLAS version of the algorithm.
 *  WARNING: Your matrix need to be diagonal dominant if you want to call this
 *  routine safely.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  @param[in] IB
 *          The block size to switch between blocked and unblocked code.
 *
 *  @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  @param[in] LDA
 *          The leading dimension of the array A.  LDA >= max(1,M).
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

int CORE_zgetrf_nopiv(int M, int N, int IB,
                      CHAMELEON_Complex64_t *A, int LDA,
                      int *INFO)
{
    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t)1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;
    int i, k, sb;
    int iinfo;

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
    if ((LDA < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    k = chameleon_min(M, N);
    for(i =0 ; i < k; i += IB) {
        sb = chameleon_min(IB, k-i);
        /*
         * Factor diagonal and subdiagonal blocks and test for exact singularity.
         */
        iinfo = CORE_zgetf2_nopiv( M-i, sb, &A[LDA*i+i], LDA );

        /*
         * Adjust info
         */
        if((*INFO == 0) && (iinfo > 0))
        	*INFO = iinfo + i;

        if (i+sb < N) {
            /*  Compute block row of U */
            CORE_ztrsm( ChamLeft, ChamLower,
                        ChamNoTrans, ChamUnit,
                        sb, N-(i+sb),
                        zone, &A[LDA*i+i],      LDA,
                              &A[LDA*(i+sb)+i], LDA);

            if ( i+sb < M )
            {
                /* Update trailing submatrix */
                CORE_zgemm(
                    ChamNoTrans, ChamNoTrans,
                    M-(i+sb), N-(i+sb), sb,
                    mzone, &A[LDA*i     + i+sb], LDA,
                           &A[LDA*(i+sb)+ i   ], LDA,
                    zone,  &A[LDA*(i+sb)+ i+sb], LDA);
            }
        }
    }
    return CHAMELEON_SUCCESS;
}
