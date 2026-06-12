/**
 *
 * @file lapack_zlange.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for lange
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-11-14
 * @precisions normal z -> s d c
 *
 */

#include "chameleon_lapack.h"
#include "lapack_api_common.h"

/* Fortran LAPACK interface */

#define CHAMELEON_lapack_zlange CHAMELEON_GLOBAL( chameleon_lapack_zlange, CHAMELEON_LAPACK_ZLANGE )
double CHAMELEON_lapack_zlange ( const char *norm, const int *m, const int *n,
                                 const CHAMELEON_Complex64_t *a, const int *lda,
                                 double *work )
{
    (void)work;
    return CHAMELEON_lapacke_zlange( CblasColMajor,
                                     *norm, *m, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlange - returns the value of the 1-norm, or the Frobenius
 *  norm, or the infinity norm, or the element of largest absolute value of a
 *  real/complex matrix A.
 *
 ******************************************************************************
 *
 * @param[in] matrix_layout Specifies whether the matrices are row or column
 *          major, it must be set to LAPACK_COL_MAJOR or CblasColMajor (102),
 *          the matrix_layout supported in Chameleon.
 *
 * @param[in] norm = 'M' or 'm': val = max(abs(Aij)), largest absolute value
 *                   of the matrix A.
 *                 = '1' or 'O' or 'o': val = norm1(A), 1-norm of the matrix A
 *                   (maximum column sum),
 *                 = 'I' or 'i': val = normI(A), infinity norm of the matrix A
 *                   (maximum row sum),
 *                 = 'F', 'f', 'E' or 'e': val = normF(A), Frobenius norm of
 *                   the matrix A (square root of sum of squares).
 *
 * @param[in] M The number of rows of the matrix A. M >= 0. When M = 0, the
 *          returned value is set to zero.
 *
 * @param[in] N The number of columns of the matrix A. N >= 0. When N = 0, the
 *          returned value is set to zero.
 *
 * @param[in] A The M-by-N matrix A.
 *
 * @param[in] LDA The leading dimension of the array A. LDA >= max(1,M).
 *
 ******************************************************************************
 *
 * @retval the norm described above.
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlange
 * @sa CHAMELEON_lapacke_clange
 * @sa CHAMELEON_lapacke_dlange
 * @sa CHAMELEON_lapacke_slange
 *
 */
double CHAMELEON_lapacke_zlange( int matrix_layout, char norm, int M, int N,
                                 const CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlange", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlange( (cham_normtype_t)chameleon_lapacktochameleon_norm(&norm),
                             M, N, (CHAMELEON_Complex64_t *)A, lda );
}
