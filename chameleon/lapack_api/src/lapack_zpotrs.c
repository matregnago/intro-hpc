/**
 *
 * @file lapack_zpotrs.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for potrs
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

#define CHAMELEON_lapack_zpotrs CHAMELEON_GLOBAL( chameleon_lapack_zpotrs, CHAMELEON_LAPACK_ZPOTRS )
void CHAMELEON_lapack_zpotrs ( const char *uplo, const int *n, const int *nrhs,
                               const CHAMELEON_Complex64_t *a, const int *lda,
                                     CHAMELEON_Complex64_t *b, const int *ldb,
                               int *info )
{
    *info = CHAMELEON_lapacke_zpotrs( CblasColMajor,
                                      *uplo, *n, *nrhs, a, *lda, b, *ldb );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zpotrs - Solves a system of linear equations A * X = B with a symmetric positive
 *  definite (or Hermitian positive definite in the complex case) matrix A using the Cholesky
 *  factorization A = U^H*U or A = L*L^H computed by CHAMELEON_zpotrf.
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
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] A
 *          The triangular factor U or L from the Cholesky factorization A = U^H*U or A = L*L^H,
 *          computed by CHAMELEON_lapacke_zpotrf.
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
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zpotrs
 * @sa CHAMELEON_lapacke_cpotrs
 * @sa CHAMELEON_lapacke_dpotrs
 * @sa CHAMELEON_lapacke_spotrs
 *
 */
int CHAMELEON_lapacke_zpotrs( int matrix_layout, char uplo, int N, int NRHS,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zpotrs", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zpotrs( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             N, NRHS, (CHAMELEON_Complex64_t *)A, lda, B, ldb );
}
