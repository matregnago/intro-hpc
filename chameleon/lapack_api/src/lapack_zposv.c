/**
 *
 * @file lapack_zposv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for posv
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

#define CHAMELEON_lapack_zposv CHAMELEON_GLOBAL( chameleon_lapack_zposv, CHAMELEON_LAPACK_ZPOSV )
void CHAMELEON_lapack_zposv ( const char *uplo, const int *n, const int *nrhs,
                              CHAMELEON_Complex64_t *a, const int *lda,
                              CHAMELEON_Complex64_t *b, const int *ldb,
                              int *info )
{
    *info = CHAMELEON_lapacke_zposv( CblasColMajor,
                                     *uplo, *n, *nrhs, a, *lda, b, *ldb );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zposv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N symmetric positive definite (or Hermitian positive definite
 *  in the complex case) matrix and X and B are N-by-NRHS matrices.
 *  The Cholesky decomposition is used to factor A as
 *
 *    \f[ A = \{_{L\times L^H, if uplo = ChamLower}^{U^H\times U, if uplo = ChamUpper} \f]
 *
 *  where U is an upper triangular matrix and  L is a lower triangular matrix.
 *  The factored form of A is then used to solve the system of equations A * X = B.
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
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If UPLO = 'U', the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the
 *          strictly lower triangular part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A
 *          contains the lower triangular part of the matrix A, and the
 *          strictly upper triangular part of A is not referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U^H*U or A = L*L^H.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, the leading minor of order i of A is not
                 positive definite, so the factorization could not be
                 completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zposv
 * @sa CHAMELEON_lapacke_cposv
 * @sa CHAMELEON_lapacke_dposv
 * @sa CHAMELEON_lapacke_sposv
 *
 */
int CHAMELEON_lapacke_zposv( int matrix_layout, char uplo, int N, int NRHS,
                             CHAMELEON_Complex64_t *A, int lda,
                             CHAMELEON_Complex64_t *B, int ldb )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zposv", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zposv( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                            N, NRHS, A, lda, B, ldb );
}
