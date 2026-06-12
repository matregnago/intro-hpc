/**
 *
 * @file core_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgemm CPU kernel
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
 * @author Guillaume Sylvand
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
 *  Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^H, \f]
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m-by-k matrix, op( B ) a k-by-n matrix and C an m-by-n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          - ChamNoTrans:   A is not transposed,
 *          - ChamTrans:     A is transposed,
 *          - ChamConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          - ChamNoTrans:   B is not transposed,
 *          - ChamTrans:     B is transposed,
 *          - ChamConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] K
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transa = ChamNoTrans,
 *          and is m otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          When transa = ChamNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transb = ChamNoTrans,
 *          and is k otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.
 *          When transb = ChamNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] LDC
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 */
void CORE_zgemm(cham_trans_t transA, cham_trans_t transB,
                int M, int N, int K,
                CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t *A, int LDA,
                                         const CHAMELEON_Complex64_t *B, int LDB,
                CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *C, int LDC)
{
#ifdef CBLAS_HAS_ZGEMM3M
  if (get_coreblas_gemm3m_enabled())
    cblas_zgemm3m(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
  else
#endif
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}


/**/








