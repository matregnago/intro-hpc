/**
 *
 * @file core_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zhemm CPU kernel
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
  *  Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *  or
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is a Hermitian matrix and B and
 *  C are m-by-n matrices.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the Hermitian matrix A appears on the
 *          left or right in the operation as follows:
 *          - ChamLeft:  \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          - ChamRight: \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the Hermitian matrix A is to be referenced as follows:
 *          - ChamLower:     Only the lower triangular part of the
 *                             Hermitian matrix A is to be referenced.
 *          - ChamUpper:     Only the upper triangular part of the
 *                             Hermitian matrix A is to be referenced.
 *
 * @param[in] M
 *          The number of rows of the matrix C. m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix C. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          A is an lda-by-ka matrix, where ka is m when side = ChamLeft,
 *          and is n otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. lda >= max(1,ka).
 *
 * @param[in] B
 *          B is an ldb-by-n matrix, where the leading m-by-n part of
 *          the array B must contain the matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          C is an ldc-by-n matrix.
 *          On exit, the array is overwritten by the m-by-n updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 */
void CORE_zhemm(cham_side_t side, cham_uplo_t uplo,
                int M, int N,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                const CHAMELEON_Complex64_t *B, int LDB,
                CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC)
{
    cblas_zhemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}
