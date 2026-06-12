/**
 *
 * @file core_zgram.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgram CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> s d c z
 *
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgram computes the Gram matrix factors of A inplace
 *
 *    \f[
 *    d_{i.}^2 = (1/n) \sum_j d_{ij}^2` and :math:`d_{..}^2 = (1/n^2) \sum_{i,j} d_{ij}^2 \\
 *    A_{i,j} = -(1/2) (d_{ij}^2 - d_{i.}^2 - d_{.j}^2 + d_{..}^2)
 *    \f]
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be computed.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in] M
 *          The number of rows of the overall matrix.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the overall matrix.  N >= 0.
 *
 * @param[in] Mt
 *          The number of rows of the tile A.  Mt >= 0.
 *
 * @param[in] Nt
 *         The number of columns of the tile A.  Nt >= 0.
 *
 * @param[in] Di
 *         The 1-by-Mt tile containing the sum of squares by rows.
 *
 * @param[in] LDDI
 *         The leading dimension of the array Di.
 *
 * @param[in] Dj
 *         The 1-by-Nt tile containing the sum of squares by columns.
 *
 * @param[in] LDDJ
 *         The leading dimension of the array Dj.
 *
 * @param[in] D
 *         The sum of squares of all the matrix.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the Gram matrix coefficients.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */

int CORE_zgram( cham_uplo_t uplo,
                int M, int N,
                int Mt, int Nt,
                const double *Di, int LDDI,
                const double *Dj, int LDDJ,
                const double *D,
                double *A, int LDA )
{
    int i, j;
    double coef = -0.5;
    double di, dj, d;

    /* Check input arguments */
    if (uplo != ChamUpper && uplo != ChamLower && uplo != ChamUpperLower ) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if (Mt < 0) {
        coreblas_error(4, "Illegal value of Mt");
        return -4;
    }
    if (Nt < 0) {
        coreblas_error(5, "Illegal value of Nt");
        return -5;
    }
    if ((LDA < chameleon_max(1,Mt)) && (Mt > 0)) {
        coreblas_error(12, "Illegal value of LDA");
        return -12;
    }

    /* Quick return */
    if ((Mt == 0) || (Nt == 0))
        return CHAMELEON_SUCCESS;

    /* overall mean of squares */
    d = ( D[0] * D[0] ) / ( (double)M * (double)N );

    for(j = 0; j < Nt; j++) {
        /* mean of squares of the column */
        dj = ( Dj[j*LDDJ] * Dj[j*LDDJ] ) / ((double)M);
        int mmin = ( uplo == ChamLower ) ? j : 0;
        int mmax = ( uplo == ChamUpper)  ? chameleon_min(j+1, Mt) : Mt;
        for(i = mmin; i < mmax; i++) {
            /* mean of squares of the row */
            di = ( Di[i*LDDI] * Di[i*LDDI] ) / ((double)N);
            /* compute Gram factor */
            A[j*LDA+i] = coef*(A[j*LDA+i]*A[j*LDA+i] - di - dj + d);
        }
    }

    return CHAMELEON_SUCCESS;
}
