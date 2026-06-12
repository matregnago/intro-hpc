/**
 *
 * @file core_ztrasm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztrasm CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_ztrasm - Computes the sums of the absolute values of elements in a same
 *  row or column in a triangular matrix.
 *  This function is an auxiliary function to triangular matrix norm computations.
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = ChamColumnwise: Computes the sum on each column
 *          = ChamRowwise:    Computes the sum on each row
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular
 *          = ChamUpper: Upper triangle of A is referenced;
 *          = ChamLower: Lower triangle of A is referenced.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: A is non unit;
 *          = ChamUnit:    A us unit.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          A is a M-by-N matrix.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension M if storev = ChamRowwise; N otherwise.
 *          On exit, contains the sums of the absolute values per column or row
 *          added to the input values.
 *
 */

void CORE_ztrasm(cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag,
                 int M, int N,
                 const CHAMELEON_Complex64_t *A, int lda, double *work)
{
    const CHAMELEON_Complex64_t *tmpA;
    int i, j, imax;
    int idiag = (diag == ChamUnit) ? 1 : 0;

    /*
     * ChamUpper / ChamColumnwise
     */
    if  (uplo == ChamUpper ) {
        M = chameleon_min(M, N);

        if (storev == ChamColumnwise) {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                imax = chameleon_min(j+1-idiag, M);

                if ( j < M )
                    work[j] += idiag;

                for (i = 0; i < imax; i++) {
                    work[j] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        /*
         * ChamUpper / ChamRowwise
         */
        else {
            if (diag == ChamUnit) {
                for (i = 0; i < M; i++) {
                    work[i] += 1.;
                }
            }
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                imax = chameleon_min(j+1-idiag, M);

                for (i = 0; i < imax; i++) {
                    work[i] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
    } else {
        N = chameleon_min(M, N);

        /*
         * ChamLower / ChamColumnwise
         */
        if (storev == ChamColumnwise) {
            for (j = 0; j < N; j++) {
                tmpA = A + j * (lda+1) + idiag;

                work[j] += idiag;
                for (i = j+idiag; i < M; i++) {
                    work[j] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        /*
         * ChamLower / ChamRowwise
         */
        else {
            if (diag == ChamUnit) {
                for (i = 0; i < N; i++) {
                    work[i] += 1.;
                }
            }
            for (j = 0; j < N; j++) {
                tmpA = A + j * (lda+1) + idiag;

                for (i = j+idiag; i < M; i++) {
                    work[i] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
    }
}
