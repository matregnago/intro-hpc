/**
 *
 * @file lapack_zpotrf.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for potrf
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

#define CHAMELEON_lapack_zpotrf CHAMELEON_GLOBAL( chameleon_lapack_zpotrf, CHAMELEON_LAPACK_ZPOTRF )
void CHAMELEON_lapack_zpotrf ( const char *uplo, const int *n,
                               CHAMELEON_Complex64_t *a, const int *lda,
                               int *info )
{
    *info = CHAMELEON_lapacke_zpotrf( CblasColMajor, *uplo, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zpotrf - Computes the Cholesky factorization of a symmetric positive definite
 *  (or Hermitian positive definite in the complex case) matrix A.
 *  The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = ChamLower}^{U^H\times U, if uplo = ChamUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = 'U' or 'u': Upper triangle of A is stored;
 *          = 'L' or 'l': Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = 'U', the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U^H*U or A = L*L^H.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zpotrf
 * @sa CHAMELEON_lapacke_cpotrf
 * @sa CHAMELEON_lapacke_dpotrf
 * @sa CHAMELEON_lapacke_spotrf
 *
 */
int CHAMELEON_lapacke_zpotrf( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zpotrf", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zpotrf( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             N, A, lda );
}
