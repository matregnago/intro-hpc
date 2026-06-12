/**
 *
 * @file core_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zher2k CPU kernel
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
 * @precisions normal z -> c
 *
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Performs one of the Hermitian rank 2k operations
 *
 *    \f[ C = \alpha A \times B^H + conjg( \alpha ) B \times A^H + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^H \times B + conjg( \alpha ) B^H \times A + \beta C, \f]
 *
 *  where alpha is a complex scalar, beta is a real scalar,
 *  C is an n-by-n Hermitian matrix, and A and B are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - ChamUpper: Upper triangle of C is stored;
 *          - ChamLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - ChamNoTrans:
 *            \f[ C = \alpha A \times B^H
 *                  + conjg( \alpha ) B \times A^H + \beta C; \f]
 *          - ChamConjTrans:
 *            \f[ C = \alpha A^H \times B
 *                  + conjg( \alpha ) B^H \times A + \beta C. \f]
 *
 * @param[in] N
 *          The order of the matrix C. n >= zero.
 *
 * @param[in] K
 *          If trans = ChamNoTrans, number of columns of the A and B matrices;
 *          if trans = ChamConjTrans, number of rows of the A and B matrices.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix.
 *          If trans = ChamNoTrans,   ka = k;
 *          if trans = ChamConjTrans, ka = n.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          If trans = ChamNoTrans,   lda >= max(1, n);
 *          if trans = ChamConjTrans, lda >= max(1, k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix.
 *          If trans = ChamNoTrans,   kb = k;
 *          if trans = ChamConjTrans, kb = n.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.
 *          If trans = ChamNoTrans,   ldb >= max(1, n);
 *          if trans = ChamConjTrans, ldb >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 */
void CORE_zher2k(cham_uplo_t uplo, cham_trans_t trans,
                 int N, int K,
                 CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                 const CHAMELEON_Complex64_t *B, int LDB,
                 double beta, CHAMELEON_Complex64_t *C, int LDC)
{
    cblas_zher2k(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        CBLAS_SADDR(alpha), A, LDA, B, LDB,
        beta, C, LDC);
}
