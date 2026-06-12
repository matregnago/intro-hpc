/**
 *
 * @file lapack_zgemm.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for gemm
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

#define CHAMELEON_blas_zgemm CHAMELEON_GLOBAL( chameleon_blas_zgemm, CHAMELEON_BLAS_ZGEMM )
void CHAMELEON_blas_zgemm ( const char *transa, const char *transb,
                            const int *m, const int *n, const int *k,
                            const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                                                                const CHAMELEON_Complex64_t *b, const int *ldb,
                            const CHAMELEON_Complex64_t *beta,        CHAMELEON_Complex64_t *c, const int *ldc )
{
    CHAMELEON_cblas_zgemm( CblasColMajor,
                           chameleon_blastocblas_trans(transa),
                           chameleon_blastocblas_trans(transb),
                           *m, *n, *k,
                           CBLAS_SADDR(*alpha), a, *lda,
                           b, *ldb,
                           CBLAS_SADDR(*beta), c, *ldc );
}

/* C CBLAS interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_cblas_zgemm - Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          Specifies whether the matrices are row or column major, it must be
 *          set to CblasColMajor, the order supported in Chameleon.
 *
 * @param[in] TransA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = CblasNoTrans:   A is not transposed;
 *          = CblasTrans:     A is transposed;
 *          = CblasConjTrans: A is conjugate transposed.
 *
 * @param[in] TransB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = CblasNoTrans:   B is not transposed;
 *          = CblasTrans:     B is transposed;
 *          = CblasConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = CblasNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = CblasNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_cblas_zgemm
 * @sa CHAMELEON_cblas_cgemm
 * @sa CHAMELEON_cblas_dgemm
 * @sa CHAMELEON_cblas_sgemm
 *
 */
void CHAMELEON_cblas_zgemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                            int M, int N, int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                               const CHAMELEON_Complex64_t *B, int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_zgemm", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
    CHAMELEON_Complex64_t betac = *(CHAMELEON_Complex64_t *)beta;
#else
    CHAMELEON_Complex64_t alphac = alpha;
    CHAMELEON_Complex64_t betac = beta;
#endif

    CHAMELEON_zgemm( (cham_trans_t)transA, (cham_trans_t)transB, M, N, K,
                     alphac, (CHAMELEON_Complex64_t *)A, lda,
                     (CHAMELEON_Complex64_t *)B, ldb,
                     betac, (CHAMELEON_Complex64_t *)C, ldc );
}
