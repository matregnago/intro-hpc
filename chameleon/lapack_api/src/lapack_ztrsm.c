/**
 *
 * @file lapack_ztrsm.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for trsm
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

#define CHAMELEON_blas_ztrsm CHAMELEON_GLOBAL( chameleon_blas_ztrsm, CHAMELEON_BLAS_ZTRSM )
void CHAMELEON_blas_ztrsm ( const char *side, const char *uplo,
                            const char *trans, const char *diag,
                            const int *m, const int *n,
                            const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                                                                      CHAMELEON_Complex64_t *b, const int *ldb )
{
    CHAMELEON_cblas_ztrsm( CblasColMajor,
                           chameleon_blastocblas_side(side),
                           chameleon_blastocblas_uplo(uplo),
                           chameleon_blastocblas_side(trans),
                           chameleon_blastocblas_diag(diag),
                           *m, *n,
                           CBLAS_SADDR(*alpha), a, *lda,
                           b, *ldb );
}

/* C CBLAS interface */

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_LAPACK_API
 *
 *  CHAMELEON_cblas_ztrsm - Solves one of the matrix equations
 *    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
 *
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
 *
 *    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
 *
 * The matrix X is overwritten on B.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          Specifies whether the matrices are row or column major, it must be
 *          set to CblasColMajor, the order supported in Chameleon.
 *
 * @param[in] side
 *          Specifies whether op(A) appears on the left or on the right of X:
 *          = CblasLeft:  op(A) * X = B
 *          = CblasRight: X * op(A) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = CblasUpper: A is an upper triangular matrix.
 *          = CblasLower: A is a lower triangular matrix.
 *
 * @param[in] trans
 *          Specifies the form of op( A ) to be used in the matrix
 *          multiplication as follows:
 *          = CblasNoTrans:   op( A ) = A.
 *          = CblasTrans:     op( A ) = A**T.
 *          = CblasConjTrans: op( A ) = A**H.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = CblasNonUnit: A is assumed to be unit triangular.
 *          = CblasUnit:    A is not assumed to be unit triangular.
 *
 * @param[in] M
 *          M specifies the number of rows of B. M must be at least zero.
 *
 * @param[in] N
 *          N specifies the number of columns of B. N must be at least zero.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha. When alpha is zero then A is not
 *          referenced and B need not be set before entry.
 *
 * @param[in] A
 *          The triangular matrix A of dimension lda-by-k,
 *          where k is M when side = CblasLeft
 *            and k is N when side = CblasRight
 *          Before entry with uplo = CblasUpper, the leading k-by-k upper
 *          triangular part of the array A must contain the upper triangular
 *          matrix and the strictly lower triangular part of A is not
 *          referenced.
 *          Before entry with uplo = CblasLower, the leading k-by-k lower
 *          triangular part of the array A must contain the lower triangular
 *          matrix and the strictly upper triangular part of A is not
 *          referenced.
 *          Note that when diag = CblasUnit, the diagonal elements of A are not
 *          referenced either, but are assumed to be unity.
 *
 * @param[in] LDA
 *          LDA specifies the first dimension of A. When side = CblasLeft then
 *          LDA must be at least max( 1, M ), when side = CblasRight then LDA
 *          must be at least max( 1, N ).
 *
 * @param[in,out] B
 *          The matrix B of dimension LDB-by-N.
 *          Before entry, the leading m by n part of the array B must contain
 *          the right-hand side matrix B, and on exit is overwritten by the
 *          solution matrix X.
 *
 * @param[in] LDB
 *          LDB specifies the first dimension of B. LDB must be at least
 *          max( 1, M ).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_cblas_ztrsm
 * @sa CHAMELEON_cblas_ctrsm
 * @sa CHAMELEON_cblas_dtrsm
 * @sa CHAMELEON_cblas_strsm
 *
 */
void CHAMELEON_cblas_ztrsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            CBLAS_TRANSPOSE trans, CBLAS_DIAG diag,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                     CHAMELEON_Complex64_t *B, int ldb )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_ztrsm", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
#else
    CHAMELEON_Complex64_t alphac = alpha;
#endif

    CHAMELEON_ztrsm( (cham_side_t)side, (cham_uplo_t)uplo,
                     (cham_trans_t)trans, (cham_diag_t)diag,
                     M, N,
                     alphac, (CHAMELEON_Complex64_t*)A, lda,
                     B, ldb );
}
