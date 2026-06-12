/**
 *
 * @file lapack_zlanhe.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for lanhe
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-11-14
 * @precisions normal z -> c
 *
 */

#include "chameleon_lapack.h"
#include "lapack_api_common.h"

/* Fortran LAPACK interface */

#define CHAMELEON_lapack_zlanhe CHAMELEON_GLOBAL( chameleon_lapack_zlanhe, CHAMELEON_LAPACK_ZLANHE )
double CHAMELEON_lapack_zlanhe ( const char *norm, const char *uplo, const int *n,
                                 const CHAMELEON_Complex64_t *a, const int *lda,
                                 double *work )
{
    (void)work;
    return CHAMELEON_lapacke_zlanhe( CblasColMajor,
                                     *norm, *uplo, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlanhe - returns the value of the one norm, or the
    Frobenius norm, or the infinity norm, or the  element of largest absolute
    value of a complex hermitian matrix A.
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
 * @param[in] uplo Specifies whether the upper or lower triangular part of the
 *                 hermitian matrix A is to be referenced.
 *                 = 'U':  Upper triangular part of A is referenced
 *                 = 'L':  Lower triangular part of A is referenced
 *
 * @param[in] N The order of the matrix A.  N >= 0.  When N = 0,
 * CHAMELEON_lapacke_zlanhe is set to zero.
 *
 * @param[in] A The hermitian matrix A, dimension (LDA,N).
 * If UPLO = 'U', the leading n by n upper triangular part of A contains the
 * upper triangular part of the matrix A, and the strictly lower triangular
 * part of A is not referenced. If UPLO = 'L', the leading n by n lower
 * triangular part of A contains the lower triangular part of the matrix A, and
 * the strictly upper triangular part of A is not referenced. Note that the
 * imaginary parts of the diagonal elements need not be set and are assumed to
 * be zero.
 *
 * @param[in] LDA The leading dimension of the array A. LDA >= max(N,1).
 *
 ******************************************************************************
 *
 * @retval the norm described above.
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlanhe
 * @sa CHAMELEON_lapacke_clanhe
 *
 */
double CHAMELEON_lapacke_zlanhe( int matrix_layout, char norm, char uplo, int N,
                                 const CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlanhe", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlanhe( (cham_normtype_t)chameleon_lapacktochameleon_norm(&norm),
                             (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             N, (CHAMELEON_Complex64_t *)A, lda );
}
