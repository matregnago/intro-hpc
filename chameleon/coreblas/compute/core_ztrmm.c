/**
 *
 * @file core_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztrmm CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
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
 *  Performs a triangular matrix-matrix multiply of the form
 *
 *          \f[B = \alpha [op(A) \times B] \f], if side = ChamLeft  or
 *          \f[B = \alpha [B \times op(A)] \f], if side = ChamRight
 *
 *  where op( X ) is one of:
 *
 *          - op(A) = A   or
 *          - op(A) = A^T or
 *          - op(A) = A^H
 *
 *  alpha is a scalar, B is an m-by-n matrix and A is a unit or non-unit, upper
 *  or lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether op( A ) appears on the left or on the right of B:
 *          - ChamLeft:  alpha*op( A )*B
 *          - ChamRight: alpha*B*op( A )
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          - ChamUpper: Upper triangle of A is stored;
 *          - ChamLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          - ChamNoTrans:   A is transposed;
 *          - ChamTrans:     A is not transposed;
 *          - ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          - ChamNonUnit: A is non-unit triangular;
 *          - ChamUnit:    A is unit triangular.
 *
 * @param[in] M
 *          The number of rows of matrix B.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of matrix B.
 *          n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A of dimension lda-by-k, where k is m when
 *          side='L' or 'l' and k is n when when side='R' or 'r'. If uplo =
 *          ChamUpper, the leading k-by-k upper triangular part of the array
 *          A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = ChamLower, the
 *          leading k-by-k lower triangular part of the array A contains the
 *          lower triangular matrix, and the strictly upper triangular part of
 *          A is not referenced. If diag = ChamUnit, the diagonal elements of
 *          A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. When side='L' or 'l',
 *          lda >= max(1,m), when side='R' or 'r' then lda >= max(1,n).
 *
 * @param[in,out] B
 *          On entry, the matrix B of dimension ldb-by-n.
 *          On exit, the result of a triangular matrix-matrix multiply
 *          ( alpha*op(A)*B ) or ( alpha*B*op(A) ).
 *
 * @param[in] LDB
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 */
void CORE_ztrmm(cham_side_t side, cham_uplo_t uplo,
                cham_trans_t transA, cham_diag_t diag,
                int M, int N,
                CHAMELEON_Complex64_t alpha,
                const CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *B, int LDB)
{
    cblas_ztrmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB);
}




