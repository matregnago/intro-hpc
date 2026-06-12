/**
 *
 * @file core_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zherk CPU kernel
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
 *  Performs one of the Hermitian rank k operations
 *
 *    \f[ C = \alpha A \times A^H + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^H \times A + \beta C, \f]
 *
 *  where alpha and beta are real scalars, C is an n-by-n Hermitian
 *  matrix, and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - ChamUpper: Upper triangle of C is stored;
 *          - ChamLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - ChamNoTrans:   \f[ C = \alpha A \times A^H + \beta C; \f]
 *          - ChamConjTrans: \f[ C = \alpha A^H \times A + \beta C. \f]
 *
 * @param[in] N
 *          The order of the matrix C. n >= 0.
 *
 * @param[in] K
 *          If trans = ChamNoTrans, number of columns of the A matrix;
 *          if trans = ChamConjTrans, number of rows of the A matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          A is an lda-by-ka matrix.
 *          If trans = ChamNoTrans,   ka = k;
 *          if trans = ChamConjTrans, ka = n.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          If trans = ChamNoTrans,   lda >= max(1, n);
 *          if trans = ChamConjTrans, lda >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          C is an ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 */
void CORE_zherk(cham_uplo_t uplo, cham_trans_t trans,
                int N, int K,
                double alpha, const CHAMELEON_Complex64_t *A, int LDA,
                double beta, CHAMELEON_Complex64_t *C, int LDC)
{
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        alpha, A, LDA,
        beta, C, LDC);
}
