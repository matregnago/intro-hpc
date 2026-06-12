/**
 *
 * @file lapack_zlauum.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lapack and lapacke api for lauum
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

#define CHAMELEON_lapack_zlauum CHAMELEON_GLOBAL( chameleon_lapack_zlauum, CHAMELEON_LAPACK_ZLAUUM )
void CHAMELEON_lapack_zlauum ( const char *uplo, const int *n,
                               CHAMELEON_Complex64_t *a, const int *lda,
                               int *info )
{
    *info = CHAMELEON_lapacke_zlauum( CblasColMajor,
                                      *uplo, *n, a, *lda );
}

/* C LAPACKE interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_lapacke_zlauum - Computes the product U * U' or L' * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.
 *
 *******************************************************************************
 *
 * @param[in] matrix_layout
 *          Specifies whether the matrices are row or column major, it must be
 *          set to LAPACK_COL_MAJOR or CblasColMajor (102), the matrix_layout
 *          supported in Chameleon.
 *
 * @param[in] uplo
 *          = 'U' or 'u':  Upper triangular
 *          = 'L' or 'l':  Lower triangular
 *
 * @param[in] N
 *          The order of the triangular factor U or L.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval =0 successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_lapacke_zlauum
 * @sa CHAMELEON_lapacke_clauum
 * @sa CHAMELEON_lapacke_dlauum
 * @sa CHAMELEON_lapacke_slauum
 *
 */
int CHAMELEON_lapacke_zlauum( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda )
{
    if ( matrix_layout != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_lapacke_zlauum", "illegal value of matrix_layout" );
        return -1;
    }

    return CHAMELEON_zlauum( (cham_uplo_t)chameleon_blastocblas_uplo(&uplo), N,
                             A, lda );
}
