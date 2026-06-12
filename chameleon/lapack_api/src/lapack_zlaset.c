/**
 *
 * @file lapack_zlaset.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for laset
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

#define CHAMELEON_lapack_zlaset CHAMELEON_GLOBAL( chameleon_lapack_zlaset, CHAMELEON_LAPACK_ZLASET )
void CHAMELEON_lapack_zlaset ( const char *uplo, const int *m, const int *n,
                               const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *beta,
                                     CHAMELEON_Complex64_t *a, const int *lda )
{
    CHAMELEON_lapacke_zlaset( CblasColMajor,
                              *uplo, *m, *n,
                              *alpha, *beta,
                              a, *lda );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlaset - The routine initializes an m-by-n matrix A to
 *  beta on the diagonal and alpha on the off-diagonals.
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A .
 *          = 'U' or 'u':  Upper triangular
 *          = 'L' or 'l':  Lower triangular
 *          Otherwise, all of the matrix A is set.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          All the off diagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlaset
 * @sa CHAMELEON_lapacke_claset
 * @sa CHAMELEON_lapacke_dlaset
 * @sa CHAMELEON_lapacke_slaset
 *
 */
int CHAMELEON_lapacke_zlaset( int matrix_layout, char uplo, int M, int N,
                              CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                              CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlaset", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlaset( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             M, N, alpha, beta, A, lda );
}
