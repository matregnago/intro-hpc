/**
 *
 * @file core_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztradd CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/lapacke.h"

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * CORE_ztradd adds two matrices together as in PBLAS pztradd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          = ChamUpperLower: A and B are general matrices.
 *          = ChamUpper: op(A) and B are upper trapezoidal matrices.
 *          = ChamLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = ChamNoTrans:   op(A) = A
 *          = ChamTrans:     op(A) = A'
 *          = ChamConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
#if defined(CHAMELEON_HAVE_WEAK)
#pragma weak CORE_ztradd = PCORE_ztradd
#define CORE_ztradd PCORE_ztradd
#define CORE_zgeadd PCORE_zgeadd
int
CORE_zgeadd(cham_trans_t trans, int M, int N,
                  CHAMELEON_Complex64_t alpha,
            const CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t beta,
                  CHAMELEON_Complex64_t *B, int LDB);
#endif
int CORE_ztradd(cham_uplo_t uplo, cham_trans_t trans, int M, int N,
                      CHAMELEON_Complex64_t  alpha,
                const CHAMELEON_Complex64_t *A, int LDA,
                      CHAMELEON_Complex64_t  beta,
                      CHAMELEON_Complex64_t *B, int LDB)
{
    int i, j, minMN;

    if (uplo == ChamUpperLower){
        int rc = CORE_zgeadd( trans, M, N, alpha, A, LDA, beta, B, LDB );
        if (rc != CHAMELEON_SUCCESS)
            return rc-1;
        else
            return rc;
    }

    if ((uplo != ChamUpper) &&
        (uplo != ChamLower))
    {
        coreblas_error(1, "illegal value of uplo");
        return -1;
    }

    if ( !isValidTrans( trans ) )
    {
        coreblas_error(2, "illegal value of trans");
        return -2;
    }

    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ( ((trans == ChamNoTrans) && (LDA < chameleon_max(1,M)) && (M > 0)) ||
         ((trans != ChamNoTrans) && (LDA < chameleon_max(1,N)) && (N > 0)) )
    {
        coreblas_error(7, "Illegal value of LDA");
        return -7;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(9, "Illegal value of LDB");
        return -9;
    }

    minMN = chameleon_min( M, N );

    if ( beta == (CHAMELEON_Complex64_t)0. ) {
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo),
                             M, N, 0., 0., B, LDB );
    }
    else if ( beta != (CHAMELEON_Complex64_t)1. ) {
        CORE_zlascal( uplo, M, N, beta, B, LDB );
    }

    /**
     * ChamLower
     */
    if (uplo == ChamLower) {
        if( trans == ChamNoTrans ) {
            for (j=0; j<minMN; j++) {
                for(i=j; i<M; i++, B++, A++) {
                    *B += alpha * (*A);
                }
                B += LDB-M+j+1;
                A += LDA-M+j+1;
            }
        }
#if defined(PRECISION_z) || defined(PRECISION_c)
        else if ( trans == ChamConjTrans ) {
            for (j=0; j<minMN; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B += alpha * conj(A[LDA*i]);
                }
                B += LDB-M+j+1;
            }
        }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
        else {
            for (j=0; j<minMN; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B += alpha * A[LDA*i];
                }
                B += LDB-M+j+1;
            }
        }
    }
    /**
     * ChamUpper
     */
    else {
        if ( trans == ChamNoTrans ) {
            for (j=0; j<N; j++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++, A++) {
                    *B += alpha * (*A);
                }
                B += LDB-mm;
                A += LDA-mm;
            }
        }
#if defined(PRECISION_z) || defined(PRECISION_c)
        else if ( trans == ChamConjTrans ) {
            for (j=0; j<N; j++, A++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B += alpha * conj(A[LDA*i]);
                }
                B += LDB-mm;
            }
        }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
        else {
            for (j=0; j<N; j++, A++) {
                int mm = chameleon_min( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B += alpha * (A[LDA*i]);
                }
                B += LDB-mm;
            }
        }
    }
    return CHAMELEON_SUCCESS;
}
