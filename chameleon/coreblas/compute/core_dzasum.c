/**
 *
 * @file core_dzasum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_dzasum CPU kernel
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
 * @defgroup CORE_CHAMELEON_Complex64_t
 * @brief Linear algebra routines exposed to users. CPU functions.
 *
 */

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_dzasum - Computes the sums of the absolute values of elements in a same
 *  row or column.
 *  This function is an auxiliary function to norm computations.
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = ChamColumnwise: Computes the sum on each column
 *          = ChamRowwise:    Computes the sum on each row
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular or general
 *          = ChamUpperLower: All matrix A is referenced;
 *          = ChamUpper: Upper triangle of A is referenced;
 *          = ChamLower: Lower triangle of A is referenced.
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
 * @param[inout] work
 *          Array of dimension M if storev = ChamRowwise; N otherwise.
 *          On exit, contains the sums of the absolute values per column or row.
 *
 */

void CORE_dzasum(cham_store_t storev, cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int lda, double *work)
{
    const CHAMELEON_Complex64_t *tmpA;
    double *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case ChamUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + cabs(*tmpA);
        }
        break;
    case ChamLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += cabs(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum;
        }
        break;
    case ChamUpperLower:
    default:
        if (storev == ChamColumnwise) {
            for (j = 0; j < N; j++) {
                /* work[j] += cblas_dzasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += cabs( A[j*lda+i] );*/
                    *tmpW += cabs( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}
