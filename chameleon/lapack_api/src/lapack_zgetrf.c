/**
 *
 * @file lapack_zgetrf.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for getrf
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-01-22
 * @precisions normal z -> s d c
 *
 */

#include "chameleon_lapack.h"
#include "lapack_api_common.h"

/* Fortran LAPACK interface */

#define CHAMELEON_lapack_zgetrf CHAMELEON_GLOBAL( chameleon_lapack_zgetrf, CHAMELEON_LAPACK_ZGETRF )
void CHAMELEON_lapack_zgetrf ( const int *m, const int *n,
                               CHAMELEON_Complex64_t *a, const int *lda,
                               int *ipiv, int *info )
{
    *info = CHAMELEON_lapacke_zgetrf( CblasColMajor, *m, *n, a, *lda, ipiv );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zgetrf - computes an LU factorization of a general M-by-N matrix A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *     A = P * L * U
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if info = i, U(i,i) is exactly zero. The factorization
 *            has been completed, but the factor U is exactly
 *            singular, and division by zero will occur if it is used
 *            to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zgetrf
 * @sa CHAMELEON_lapacke_cgetrf
 * @sa CHAMELEON_lapacke_dgetrf
 * @sa CHAMELEON_lapacke_sgetrf
 *
 */
int CHAMELEON_lapacke_zgetrf( int matrix_layout, int M, int N,
                              CHAMELEON_Complex64_t *A, int lda,
                              int *IPIV )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zgetrf", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zgetrf( M, N, A, lda, IPIV );
}
