/**
 *
 * @file lapack_ztrmm.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas and cblas api for trmm
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

#define CHAMELEON_blas_ztrmm CHAMELEON_GLOBAL( chameleon_blas_ztrmm, CHAMELEON_BLAS_ZTRMM )
void CHAMELEON_blas_ztrmm ( const char *side, const char *uplo,
                            const char *trans, const char *diag,
                            const int *m, const int *n,
                            const CHAMELEON_Complex64_t *alpha, const CHAMELEON_Complex64_t *a, const int *lda,
                                                                      CHAMELEON_Complex64_t *b, const int *ldb )
{
    CHAMELEON_cblas_ztrmm( CblasColMajor,
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
 *  CHAMELEON_cblas_ztrmm - Computes B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *******************************************************************************
 *
 * @param[in] order
 *          Specifies whether the matrices are row or column major, it must be
 *          set to CblasColMajor, the order supported in Chameleon.
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = CblasLeft:  A*X = B
 *          = CblasRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = CblasUpper: Upper triangle of A is stored;
 *          = CblasLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = CblasNoTrans:   A is transposed;
 *          = CblasTrans:     A is not transposed;
 *          = CblasConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = CblasNonUnit: A is non unit;
 *          = CblasUnit:    A us unit.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B if side = CblasLeft or
 *          the number of rows of the matrix B if side = CblasRight. NRHS >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = CblasUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = CblasLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = CblasUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X if side = CblasLeft
 *          or the NRHS-by-N matrix X if side = CblasRight.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_cblas_ztrmm
 * @sa CHAMELEON_cblas_ctrmm
 * @sa CHAMELEON_cblas_dtrmm
 * @sa CHAMELEON_cblas_strmm
 *
 */
void CHAMELEON_cblas_ztrmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            CBLAS_TRANSPOSE trans, CBLAS_DIAG diag,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                     CHAMELEON_Complex64_t *B, int ldb )
{
    if ( order != CblasColMajor ){
        fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", "CHAMELEON_cblas_ztrmm", "illegal value of order" );
        return;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    CHAMELEON_Complex64_t alphac = *(CHAMELEON_Complex64_t *)alpha;
#else
    CHAMELEON_Complex64_t alphac = alpha;
#endif

    CHAMELEON_ztrmm( (cham_side_t)side, (cham_uplo_t)uplo,
                     (cham_trans_t)trans, (cham_diag_t)diag,
                     M, N,
                     alphac, (CHAMELEON_Complex64_t *)A, lda,
                     B, ldb );
}
