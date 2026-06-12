/**
 *
 * @file core_ztplqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztplqt CPU kernel
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
 * CORE_ztplqt computes a blocked LQ factorization of a complex
 * "triangular-pentagonal" matrix C, which is composed of a
 * triangular block A and pentagonal block B, using the compact
 * WY representation for Q.
 *
 *  C = | A B | = L * Q
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile B, and the order of the
 *         triangular matrix A. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix B.
 *
 * @param[in] L
 *          The number of rows of the lower trapezoidal part of B.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in] IB
 *         The inner-blocking size.  M >= IB >= 0.
 *
 * @param[in,out] A
 *          On entry, the lower triangular M-by-M matrix A.
 *          On exit, the elements on and below the diagonal of the array
 *          contain the lower triangular matrix L.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] B
 *          On entry, the pentagonal M-by-N matrix B. The first N-L columns
 *          are rectangular, and the last L columns are lower trapezoidal.
 *          On exit, B contains the pentagonal matrix V. See Further Details.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.  LDB >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         The lower triangular block reflectors stored in compact form
 *         as a sequence of upper triangular blocks. See Further Details.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *          WORK is COMPLEX*16 array, dimension ((IB+1)*M)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int CORE_ztplqt( int M, int N, int L, int IB,
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
    if ((LDA < chameleon_max(1,M)) && (M > 0)) {
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
        //LAPACKE_ztpmlqt_work( LAPACK_COL_MAJOR, M, N, K, L, IB, V, LDV, T, LDT, A, LDA, B, LDB, WORK );
        coreblas_error( 6, "Illegal value of L (only 0 or min(M,N) handled for now)");
        return -6;
    }
#endif /*!defined(NDEBUG)*/

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return CHAMELEON_SUCCESS;

    if ( L == 0 ) {
        CORE_ztslqt( M, N, IB, A, LDA, B, LDB, T, LDT, WORK, WORK+M );
    }
    else /* if (L == M) */ {
        CORE_zttlqt( M, N, IB, A, LDA, B, LDB, T, LDT, WORK, WORK+M );
    }
    /* else { */
    /*     //LAPACKE_ztplqt_work( LAPACK_COL_MAJOR, M, N, L, IB, A, LDA, B, LDB, T, LDT, WORK ); */
    /*     coreblas_error( 3, "Illegal value of L (only 0 or M handled for now)"); */
    /*     return -3; */
    /* } */
    return CHAMELEON_SUCCESS;
}
