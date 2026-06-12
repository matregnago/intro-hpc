/**
 *
 * @file core_zgesum.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgesum CPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 * @brief Subcase storev == ChamColumnwise of CORE_zgesum()
 */
static inline int
CORE_zgesum_col( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *sum )
{
    int i, j;
    const CHAMELEON_Complex64_t *ptr = A;
    for(j=0; j<N; j++) {
        for(i=0; i<M; i++) {
            sum[j] += *ptr;
            ptr++;
        }
        ptr += LDA - M;
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamRowwise of CORE_zgesum()
 */
static inline int
CORE_zgesum_row( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *sum )
{
    int i, j;
    const CHAMELEON_Complex64_t *ptr = A;
    for(j=0; j<N; j++) {
        for(i=0; i<M; i++) {
            sum[i] += *ptr;
            ptr++;
        }
        ptr += LDA - M;
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamEltwise of CORE_zgesum()
 */
static inline int
CORE_zgesum_elt( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *sum )
{
    int i, j;
    const CHAMELEON_Complex64_t *ptr = A;
    for(j=0; j<N; j++) {
        for(i=0; i<M; i++) {
            sum[0] += *ptr;
            ptr++;
        }
        ptr += LDA - M;
    }
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgesum returns the sum of a matrix values over a chosen axis:
 *  column-wise, row-wise, element-wise.
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = ChamColumnwise: Computes the sum of values on each column
 *          = ChamRowwise:    Computes the sum of values on each row
 *          = ChamEltwise:    Computes the sum of values on all the matrix
 *
 *  @param[in] M
 *          The number of rows in the tile A.
 *
 *  @param[in] N
 *          The number of columns in the tile A.
 *
 *  @param[in] A
 *          The M-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,M).
 *
 *  @param[out] sum resulting sums
 *          Dimension:  (1,K)
 *          if storev == ChamColumnwise, K = N
 *          if storev == ChamRowwise,    K = M
 *          if storev == ChamEltwise,    K = 1
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval -k, the k-th argument had an illegal value
 *
 */
int CORE_zgesum( cham_store_t storev, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *sum )
{
    int i;
    int K;

    /* initialize pairs scale, sumsquare if not already done */
    if ( storev == ChamColumnwise ) {
        K = N;
    } else if ( storev == ChamRowwise ) {
        K = M;
    } else {
        K = 1;
    }
    for (i=0; i<K; i+=1) {
        sum[i] = 0.;
    }

    if (storev == ChamColumnwise) {
        CORE_zgesum_col( M, N, A, LDA, sum );
    }
    else if (storev == ChamRowwise) {
        CORE_zgesum_row( M, N, A, LDA, sum );
    }
    else {
        CORE_zgesum_elt( M, N, A, LDA, sum );
    }
    return CHAMELEON_SUCCESS;
}
