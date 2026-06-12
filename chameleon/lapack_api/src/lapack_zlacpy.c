/**
 *
 * @file lapack_zlacpy.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for lacpy
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

#define CHAMELEON_lapack_zlacpy CHAMELEON_GLOBAL( chameleon_lapack_zlacpy, CHAMELEON_LAPACK_ZLACPY )
void CHAMELEON_lapack_zlacpy ( const char *uplo, const int *m, const int *n,
                               const CHAMELEON_Complex64_t *a, const int *lda,
                                     CHAMELEON_Complex64_t *b, const int *ldb )
{
    CHAMELEON_lapacke_zlacpy( CblasColMajor,
                              *uplo, *m, *n,
                              a, *lda,
                              b, *ldb );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlacpy - copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *          = 'U' or 'u':  Upper triangular
 *          = 'L' or 'l':  Lower triangular
 *          Otherwise, all of the matrix A is copied.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = 'U', only the upper triangle or
 *          trapezoid is accessed; if uplo = 'L', only the lower triangle or
 *          trapezoid is accessed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by uplo.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlacpy
 * @sa CHAMELEON_lapacke_clacpy
 * @sa CHAMELEON_lapacke_dlacpy
 * @sa CHAMELEON_lapacke_slacpy
 *
 */
int CHAMELEON_lapacke_zlacpy( int matrix_layout, char uplo, int M, int N,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlacpy", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlacpy( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             M, N, (CHAMELEON_Complex64_t *)A, lda, B, ldb );
}
