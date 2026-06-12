/**
 *
 * @file core_ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztrsm CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
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
 *  Solves one of the matrix equations
 *
 *    \f[ op( A )\times X  = \alpha B, \f] or
 *    \f[ X \times op( A ) = \alpha B, \f]
 *
 *  where op( A ) is one of:
 *    \f[ op( A ) = A,   \f]
 *    \f[ op( A ) = A^T, \f]
 *    \f[ op( A ) = A^H, \f]
 *
 *  alpha is a scalar, X and B are m-by-n matrices, and
 *  A is a unit or non-unit, upper or lower triangular matrix.
 *  The matrix X overwrites B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          - ChamLeft:  op(A)*X = B,
 *          - ChamRight: X*op(A) = B.
 *
 * @param[in] uplo
 *          - ChamUpper: A is upper triangular,
 *          - ChamLower: A is lower triangular.
 *
 * @param[in] transA
 *          - ChamNoTrans:   A is not transposed,
 *          - ChamTrans:     A is transposed,
 *          - ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          - ChamNonUnit: A has non-unit diagonal,
 *          - ChamUnit:    A has unit diagonal.
 *
 * @param[in] M
 *          The number of rows of the matrix B. m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix B. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The lda-by-ka triangular matrix,
 *          where ka = m if side = ChamLeft,
 *            and ka = n if side = ChamRight.
 *          If uplo = ChamUpper, the leading k-by-k upper triangular part
 *          of the array A contains the upper triangular matrix, and the
 *          strictly lower triangular part of A is not referenced.
 *          If uplo = ChamLower, the leading k-by-k lower triangular part
 *          of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced.
 *          If diag = ChamUnit, the diagonal elements of A are also not
 *          referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[in,out] B
 *          On entry, the ldb-by-n right hand side matrix B.
 *          On exit, if return value = 0, the ldb-by-n solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 */
void CORE_ztrsm(cham_side_t side, cham_uplo_t uplo,
                cham_trans_t transA, cham_diag_t diag,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *B, int LDB)
{
    cblas_ztrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB);
}


