/**
 *
 * @file core_ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztpqrt CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * CORE_ztpqrt computes a blocked QR factorization of a complex
 * "triangular-pentagonal" matrix C, which is composed of a
 * triangular block A and pentagonal block B, using the compact
 * WY representation for Q.
 *
 *  C = | A | = Q * R
 *      | B |
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile B. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix B, and the order of the matrix
 *          A. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  N >= IB >= 0.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of B.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in,out] A
 *          On entry, the upper triangular N-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array
 *          contain the upper triangular matrix R.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
 *          are rectangular, and the last L rows are upper trapezoidal.
 *          On exit, B contains the pentagonal matrix V.  See Further Details.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.  LDB >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *          WORK is COMPLEX*16 array, dimension ((IB+1)*N)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int CORE_ztpqrt( int M, int N, int L, int IB,
                 CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB,
                 CHAMELEON_Complex64_t *T, int LDT,
                 CHAMELEON_Complex64_t *WORK )
{
#if !defined(NDEBUG)
     /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if( (L < 0) || ((L > chameleon_min(M, N)) && (chameleon_min(M,N) > 0))) {
        coreblas_error(3, "Illegal value of L");
        return -3;
    }
    if (IB < 0) {
        coreblas_error(4, "Illegal value of IB");
        return -4;
    }
    if ((LDA < chameleon_max(1,N)) && (N > 0)) {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ((LDB < chameleon_max(1,M)) && (M > 0)) {
        coreblas_error(6, "Illegal value of LDB");
        return -8;
    }
    if ((LDT < chameleon_max(1,IB)) && (IB > 0)) {
        coreblas_error(6, "Illegal value of LDT");
        return -10;
    }
    if ((L != 0) && (L != chameleon_min(M, N))) {
        //LAPACKE_ztpmqrt_work( LAPACK_COL_MAJOR, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK );
        coreblas_error( 6, "Illegal value of L (only 0 or min(M,N) handled for now)");
        return -6;
    }
#endif /*!defined(NDEBUG)*/

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    if ( L == 0 ) {
        CORE_ztsqrt( M, N, IB, A, LDA, B, LDB, T, LDT, WORK, WORK+N );
    }
    else /* if (L == M) */ {
        CORE_zttqrt( M, N, IB, A, LDA, B, LDB, T, LDT, WORK, WORK+N );
    }
    /* else { */
    /*     //LAPACKE_ztpqrt_work( LAPACK_COL_MAJOR, M, N, L, IB, A, LDA, B, LDB, T, LDT, WORK ); */
    /*     coreblas_error( 3, "Illegal value of L (only 0 or M handled for now)"); */
    /*     return -3; */
    /* } */
    return CHAMELEON_SUCCESS;
}
