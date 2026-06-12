/**
 *
 * @file core_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgeadd CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Florent Pruvost
 * @author Nicolas Brieuc
 * @date 2024-11-26
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
 *  CORE_zgeadd adds to matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = ChamNoTrans:   op(A) = A
 *          = ChamTrans:     op(A) = A'
 *          = ChamConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = ChamNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = ChamNoTrans, and k=N otherwise.
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
#pragma weak CORE_zgeadd = PCORE_zgeadd
#define CORE_zgeadd PCORE_zgeadd
#endif
int CORE_zgeadd(cham_trans_t trans, int M, int N,
                      CHAMELEON_Complex64_t  alpha,
                const CHAMELEON_Complex64_t *A, int LDA,
                      CHAMELEON_Complex64_t  beta,
                      CHAMELEON_Complex64_t *B, int LDB)
{
    int i, j;

    if ( !isValidTrans( trans ) )
    {
        coreblas_error(1, "illegal value of trans");
        return -1;
    }

    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if ( ((trans == ChamNoTrans) && (LDA < chameleon_max(1,M)) && (M > 0)) ||
         ((trans != ChamNoTrans) && (LDA < chameleon_max(1,N)) && (N > 0)) )
    {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    if ( beta == 0. ) {
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A',
                             M, N, 0., 0., B, LDB );
    }
    else if ( beta != 1. ) {
        CORE_zlascal( ChamUpperLower, M, N, beta, B, LDB );
    }

    if( trans == ChamNoTrans ) {
        for (j=0; j<N; j++) {
            for(i=0; i<M; i++, B++, A++) {
                *B += alpha * (*A);
            }
            A += LDA-M;
            B += LDB-M;
        }
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    else if ( trans == ChamConjTrans ) {
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B += alpha * conj(A[LDA*i]);
            }
            B += LDB-M;
        }
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
    else {
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B += alpha * A[LDA*i];
            }
            B += LDB-M;
        }
    }

    return CHAMELEON_SUCCESS;
}
