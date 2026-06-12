/**
 *
 * @file lapack_zsymm.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for symm
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

#define CHAMELEON_blas_zsymm CHAMELEON_GLOBAL( chameleon_blas_zsymm, CHAMELEON_BLAS_ZSYMM )
void CHAMELEON_blas_zsymm ( const char *side, const char *uplo,
                            const int *m, const int *n,
                            const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                                                                const CHAMELEON_Complex64_t *b, const int *ldb,
                            const CHAMELEON_Complex64_t *beta,  CHAMELEON_Complex64_t *c, const int *ldc )
{
    CHAMELEON_cblas_zsymm( CblasColMajor,
                           chameleon_blastocblas_side(side),
                           chameleon_blastocblas_uplo(uplo),
                           *m, *n,
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
 *  CHAMELEON_cblas_zsymm - Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *
 *  or
 *
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is a symmetric matrix and  B and
 *  C are m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          Specifies whether the matrices are row or column major, it must be
 *          set to CblasColMajor, the order supported in Chameleon.
 *
 * @param[in] side
 *          Specifies whether the symmetric matrix A appears on the
 *          left or right in the operation as follows:
 *          = CblasLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = CblasRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = CblasLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = CblasUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 * @param[in] M
 *          Specifies the number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          Specifies the number of columns of the matrix C. N >= 0.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = CblasLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,ka).
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_cblas_zsymm
 * @sa CHAMELEON_cblas_csymm
 * @sa CHAMELEON_cblas_dsymm
 * @sa CHAMELEON_cblas_ssymm
 *
 */
void CHAMELEON_cblas_zsymm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                               const CHAMELEON_Complex64_t *B, int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_zsymm", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
    CHAMELEON_Complex64_t betac = *(CHAMELEON_Complex64_t *)beta;
#else
    CHAMELEON_Complex64_t alphac = alpha;
    CHAMELEON_Complex64_t betac = beta;
#endif

    CHAMELEON_zsymm( (cham_side_t)side, (cham_uplo_t)uplo, M, N,
                     alphac, (CHAMELEON_Complex64_t *)A, lda,
                     (CHAMELEON_Complex64_t *)B, ldb,
                     betac, C, ldc );
}
