/**
 *
 * @file lapack_zlantr.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for lantr
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

#define CHAMELEON_lapack_zlantr CHAMELEON_GLOBAL( chameleon_lapack_zlantr, CHAMELEON_LAPACK_ZLANTR )
double CHAMELEON_lapack_zlantr ( const char *norm, const char *uplo, const char *diag,
                                 const int *m, const int *n,
                                 const CHAMELEON_Complex64_t *a, const int *lda,
                                 double *work )
{
    (void)work;
    return CHAMELEON_lapacke_zlantr( CblasColMajor,
                                     *norm, *uplo, *diag, *m, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlantr - returns the value of the one norm,  or the
    Frobenius norm, or the  infinity norm,  or the  element of  largest absolute
    value  of a trapezoidal or triangular matrix A.
 *
 ******************************************************************************
 *
 * @param[in] matrix_layout Specifies whether the matrices are row or column
 *          major, it must be set to LAPACK_COL_MAJOR or CblasColMajor (102),
 *          the matrix_layout supported in Chameleon.
 *
 * @param[in] norm = max(abs(A(i,j))), NORM = 'M' or 'm'
 *                 = norm1(A),         NORM = '1', 'O' or 'o'
 *                 = normI(A),         NORM = 'I' or 'i'
 *                 = normF(A),         NORM = 'F', 'f', 'E' or 'e'
 *
 * where norm1 denotes the one norm of a matrix (maximum column sum), normI
 * denotes the infinity norm of a matrix (maximum row sum) and normF denotes
 * the Frobenius norm of a matrix (square root of sum of squares). Note that
 * max(abs(A(i,j))) is not a consistent matrix norm.
 *
 * @param[in] uplo Specifies whether the matrix A is upper or lower trapezoidal.
 *                 = 'U':  Upper trapezoidal
 *                 = 'L':  Lower trapezoidal
 *                 Note that A is triangular instead of trapezoidal if M = N.
 *
 * @param[in] diag Specifies whether or not the matrix A has unit diagonal.
 *                 = 'N':  Non-unit diagonal
 *                 = 'U':  Unit diagonal
 *
 * @param[in] M The number of rows of the matrix A. M >= 0. When M = 0, the
 *          returned value is set to zero.
 *
 * @param[in] N The number of columns of the matrix A. N >= 0. When N = 0, the
 *          returned value is set to zero.
 *
 * @param[in] A The trapezoidal matrix A (A is triangular if M = N),
 * dimension (LDA,N).
 * If UPLO = 'U', the leading m by n upper trapezoidal part of
 * the array A contains the upper trapezoidal matrix, and the
 * strictly lower triangular part of A is not referenced.
 * If UPLO = 'L', the leading m by n lower trapezoidal part of
 * the array A contains the lower trapezoidal matrix, and the
 * strictly upper triangular part of A is not referenced.  Note
 * that when DIAG = 'U', the diagonal elements of A are not
 * referenced and are assumed to be one.
 *
 * @param[in] LDA The leading dimension of the array A. LDA >= max(M,1).
 *
 ******************************************************************************
 *
 * @retval the norm described above.
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlantr
 * @sa CHAMELEON_lapacke_clantr
 * @sa CHAMELEON_lapacke_dlantr
 * @sa CHAMELEON_lapacke_slantr
 *
 */
double CHAMELEON_lapacke_zlantr( int matrix_layout, char norm, char uplo, char diag,
                                 int M, int N, const CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlantr", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlantr( (cham_normtype_t)chameleon_lapacktochameleon_norm(&norm),
                             (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             (cham_diag_t)chameleon_blastocblas_diag(&diag),
                             M, N, (CHAMELEON_Complex64_t *)A, lda );
}
