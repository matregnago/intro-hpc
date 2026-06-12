/**
 *
 * @file lapack_zsyrk.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for syrk
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

/* Fortran BLAS interface */

#define CHAMELEON_blas_zsyrk CHAMELEON_GLOBAL( chameleon_blas_zsyrk, CHAMELEON_BLAS_ZSYRK )
void CHAMELEON_blas_zsyrk ( const char *uplo, const char *trans,
                            const int *n, const int *k,
                            const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                            const CHAMELEON_Complex64_t *beta,        CHAMELEON_Complex64_t *c, const int *ldc )
{
    CHAMELEON_cblas_zsyrk( CblasColMajor,
                           chameleon_blastocblas_uplo(uplo),
                           chameleon_blastocblas_trans(trans),
                           *n, *k,
                           CBLAS_SADDR(*alpha), a, *lda,
                           CBLAS_SADDR(*beta), c, *ldc );
}

/* C CBLAS interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_cblas_zsyrk - Performs one of the symmetric rank k operations
 *
 *    \f[ C = \alpha [ A \times A' ] + \beta C \f]
 *    or
 *    \f[ C = \alpha [ A' \times A ] + \beta C \f]
 *
 *  where alpha and beta are real scalars, C is an n-by-n symmetric
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
 *          Specifies the operation to be performed as follows:
 *          = CblasNoTrans: \f[ C = \alpha [ A \times A' ] + \beta C \f]
 *          = CblasTrans: \f[ C = \alpha [ A' \times A ] + \beta C \f]
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix A.
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
 * @sa CHAMELEON_cblas_zsyrk
 * @sa CHAMELEON_cblas_csyrk
 * @sa CHAMELEON_cblas_dsyrk
 * @sa CHAMELEON_cblas_ssyrk
 *
 */
void CHAMELEON_cblas_zsyrk( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                            int N, int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_zsyrk", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
    CHAMELEON_Complex64_t betac = *(CHAMELEON_Complex64_t *)beta;
#else
    CHAMELEON_Complex64_t alphac = alpha;
    CHAMELEON_Complex64_t betac = beta;
#endif

    CHAMELEON_zsyrk( (cham_uplo_t)uplo, (cham_trans_t)trans, N, K,
                     alphac, (CHAMELEON_Complex64_t *)A, lda,
                     betac, C, ldc );
}
