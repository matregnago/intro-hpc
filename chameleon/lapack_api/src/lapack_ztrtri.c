/**
 *
 * @file lapack_ztrtri.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for trtri
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

#define CHAMELEON_lapack_ztrtri CHAMELEON_GLOBAL( chameleon_lapack_ztrtri, CHAMELEON_LAPACK_ZTRTRI )
void CHAMELEON_lapack_ztrtri ( const char *uplo, const char *diag, const int *n,
                               CHAMELEON_Complex64_t *a, const int *lda,
                               int *info )
{
    *info = CHAMELEON_lapacke_ztrtri( CblasColMajor,
                                      *uplo, *diag, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_ztrtri - Computes the inverse of a complex upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] uplo
 *          = 'U' or 'u': Upper triangle of A is stored;
 *          = 'L' or 'l': Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = 'N', then A is not a unit triangular matrix;
 *          = 'U', A is unit triangular: diagonal elements of A are assumed to
 *            be 1 and not referenced in the array a.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If UPLO =
 *          'L', the leading N-by-N lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If DIAG =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix, in the same storage format.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, A(i,i) is exactly zero.  The triangular
 *               matrix is singular and its inverse can not be computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_ztrtri
 * @sa CHAMELEON_lapacke_ctrtri
 * @sa CHAMELEON_lapacke_dtrtri
 * @sa CHAMELEON_lapacke_strtri
 *
 */
int CHAMELEON_lapacke_ztrtri( int matrix_layout, char uplo, char diag, int N,
                              CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_ztrtri", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_ztrtri( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo),
                             (cham_diag_t)chameleon_blastocblas_diag(&diag),
                             N, A, lda );
}
