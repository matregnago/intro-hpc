/**
 *
 * @file lapack_zher2k.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for her2k
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

#define CHAMELEON_blas_zher2k CHAMELEON_GLOBAL( chameleon_blas_zher2k, CHAMELEON_BLAS_ZHER2K )
void CHAMELEON_blas_zher2k ( const char *uplo, const char *trans,
                             const int *n, const int *k,
                             const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                                                                 const CHAMELEON_Complex64_t *b, const int *ldb,
                             const double *beta,                       CHAMELEON_Complex64_t *c, const int *ldc )
{
    CHAMELEON_cblas_zher2k( CblasColMajor,
                            chameleon_blastocblas_uplo(uplo),
                            chameleon_blastocblas_trans(trans),
                            *n, *k,
                            CBLAS_SADDR(*alpha), a, *lda,
                            b, *ldb,
                            *beta, c, *ldc );
}

/* C CBLAS interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_cblas_zher2k - Performs one of the hermitian rank 2k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( B )' )] + conjg( \alpha ) [ op( B ) \times conjg( op( A )' )] + \beta C \f],
 *    or
 *    \f[ C = \alpha [ conjg( op( A )' ) \times op( B ) ] + conjg( \alpha ) [ conjg( op( B )' ) \times op( A ) ] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A and B are an n-by-k matrices the first case and k-by-n
 *  matrices in the second case.
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
 *          = CblasNoTrans:   \f[ C = \alpha [ op( A ) \times conjg( op( B )' )] + conjg( \alpha ) [ op( B ) \times conjg( op( A )' )] + \beta C \f]
 *          = CblasConjTrans: \f[ C = \alpha [ conjg( op( A )' ) \times op( B ) ] + conjg( \alpha ) [ conjg( op( B )' ) \times op( A ) ] + \beta C \f]
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the A and B matrices with trans = CblasNoTrans.
 *          K specifies the number of rows of the A and B matrices with trans = CblasTrans.
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
 *          max( 1, N ), otherwise LDA must be at least max( 1, K ).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is K when trans = CblasNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB must be at least
 *          max( 1, N ), otherwise LDB must be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
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
 * @sa CHAMELEON_cblas_zher2k
 * @sa CHAMELEON_cblas_cher2k
 * @sa CHAMELEON_cblas_dher2k
 * @sa CHAMELEON_cblas_sher2k
 *
 */
void CHAMELEON_cblas_zher2k( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                             int N, int K,
                             const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                const CHAMELEON_Complex64_t *B, int ldb,
                             double beta,             CHAMELEON_Complex64_t *C, int ldc )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_zher2k", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
#else
    CHAMELEON_Complex64_t alphac = alpha;
#endif

    CHAMELEON_zher2k( (cham_uplo_t)uplo, (cham_trans_t)trans, N, K,
                      alphac, (CHAMELEON_Complex64_t *)A, lda,
                      (CHAMELEON_Complex64_t *)B, ldb,
                      beta, (CHAMELEON_Complex64_t *)C, ldc );
}
