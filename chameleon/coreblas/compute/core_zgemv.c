/**
 *
 * @file core_zgemv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgemv CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Performs one of the matrix-vector operations
 *
 *    \f[ y = \alpha [op( A ) \times x] + \beta y \f],
 *
 *  where op( A ) is one of:
 *    \f[ op( A ) = A,   \f]
 *    \f[ op( A ) = A^T, \f]
 *    \f[ op( A ) = A^H, \f]
 *
 *  alpha and beta are scalars, op(A) an m-by-n matrix, and x and y are two vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          - ChamNoTrans:   A is not transposed,
 *          - ChamTrans:     A is transposed,
 *          - ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-n matrix, where only the m-by-n leading entries are references.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] x
 *          X is COMPLEX*16 array, dimension at least
 *          ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
 *          and at least
 *          ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
 *          Before entry, the incremented array X must contain the
 *          vector x.
 *
 * @param[in] incX
 *          On entry, INCX specifies the increment for the elements of
 *          X. INCX must not be zero.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] y
 *          Y is COMPLEX*16 array, dimension at least
 *          ( 1 + ( n - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
 *          and at least
 *          ( 1 + ( m - 1 )*abs( INCY ) ) otherwise.
 *          Before entry, the incremented array Y must contain the vector y with
 *          beta != 0. On exit, Y is overwritten by the updated vector.
 *
 * @param[in] incY
 *          On entry, INCY specifies the increment for the elements of
 *          Y. INCY must not be zero.
 *
 */
void CORE_zgemv( cham_trans_t trans, int M, int N,
                 CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                              const CHAMELEON_Complex64_t *x, int incX,
                 CHAMELEON_Complex64_t beta,        CHAMELEON_Complex64_t *y, int incY )
{
    cblas_zgemv(
        CblasColMajor, (CBLAS_TRANSPOSE)trans, M, N,
        CBLAS_SADDR(alpha), A, LDA,
                            x, incX,
        CBLAS_SADDR(beta),  y, incY );
}
