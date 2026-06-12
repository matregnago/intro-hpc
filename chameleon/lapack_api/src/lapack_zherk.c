/**
 *
 * @file lapack_zherk.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for herk
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

/* Fortran BLAS interface */

#define CHAMELEON_blas_zherk CHAMELEON_GLOBAL( chameleon_blas_zherk, CHAMELEON_BLAS_ZHERK )
void CHAMELEON_blas_zherk ( const char *uplo, const char *trans,
                            const int *n, const int *k,
                            const double *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                            const double *beta,        CHAMELEON_Complex64_t *c, const int *ldc )
{
    CHAMELEON_cblas_zherk( CblasColMajor,
                           chameleon_blastocblas_uplo(uplo),
                           chameleon_blastocblas_trans(trans),
                           *n, *k,
                           *alpha, a, *lda,
                           *beta, c, *ldc );
}

/* C CBLAS interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_cblas_zherk - Performs one of the hermitian rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          Specifies whether the matrices are row or column major, it must be
 *          set to CblasColMajor, the order supported in Chameleon.
 *
 * @param[in] uplo
 *          = CblasUpper: Upper triangle of C is stored;
 *          = CblasLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = CblasNoTrans:   A is not transposed;
 *          = CblasConjTrans: A is conjugate transposed.
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = CblasNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ) if trans == CblasNoTrans, otherwise LDA must
 *          be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_cblas_zherk
 * @sa CHAMELEON_cblas_cherk
 * @sa CHAMELEON_cblas_dherk
 * @sa CHAMELEON_cblas_sherk
 *
 */
void CHAMELEON_cblas_zherk( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                            int N, int K,
                            double alpha, const CHAMELEON_Complex64_t *A, int lda,
                            double beta,        CHAMELEON_Complex64_t *C, int ldc )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_zherk", "illegal value of order" );
        return;
    }

    CHAMELEON_zherk( (cham_uplo_t)uplo, (cham_trans_t)trans, N, K,
                     alpha, (CHAMELEON_Complex64_t *)A, lda,
                     beta, (CHAMELEON_Complex64_t *)C, ldc );
}
