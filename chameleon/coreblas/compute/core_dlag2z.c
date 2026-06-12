/**
 *
 * @file core_dlag2z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_dlag2z CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Converts a m-by-n matrix A from double real to double complex.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *          n >= 0.
 *
 * @param[in] A
 *          The lda-by-n matrix in double complex precision to convert.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          On exit, the converted LDB-by-n matrix in single complex precision.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix As.
 *          ldas >= max(1,m).
 *
 */
void
CORE_dlag2z( cham_uplo_t uplo, int M, int N,
             const double          *A, int lda,
             CHAMELEON_Complex64_t *B, int ldb )
{
    const double          *Aptr;
    CHAMELEON_Complex64_t *Bptr;
    int i, j;

    if ( (uplo != ChamUpperLower) &&
         (uplo != ChamUpper)      &&
         (uplo != ChamLower))
    {
        coreblas_error(1, "illegal value of uplo");
        return;
    }

    if (M < 0) {
        coreblas_error(2, "Illegal value of m");
        return;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of n");
        return;
    }
    if ( (lda < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(5, "Illegal value of lda");
        return;
    }
    if ( (ldb < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(7, "Illegal value of ldb");
        return;
    }

    for(j=0; j<N; j++) {
        int mmin = ( uplo == ChamLower ) ? j                     : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(j+1, M) : M;

        Aptr = A + lda * j + mmin;
        Bptr = B + ldb * j + mmin;

        for(i=mmin; i<mmax; i++, Aptr++, Bptr++) {
            *Bptr = (CHAMELEON_Complex64_t)(*Aptr);
        }
    }
}
