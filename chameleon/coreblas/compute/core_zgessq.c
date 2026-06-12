/**
 *
 * @file core_zgessq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgessq CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "coreblas/lapacke.h"
#include "coreblas/sumsq_update.h"
#include "coreblas.h"

/**
 * @brief Subcase storev == ChamColumnwise of CORE_zgessq()
 */
static inline int
CORE_zgessq_col( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq )
{
    int i, j;
    double tmp;
    double *ptr, *tmpScale, *tmpSumsq;

    for(j=0; j<N; j++) {
        ptr = (double*) ( A + j * LDA );
        tmpScale = sclssq+2*j;
        tmpSumsq = sclssq+2*j+1;
        for(i=0; i<M; i++, ptr++) {
            tmp = fabs(*ptr);
            sumsq_update( 1., tmpScale, tmpSumsq, &tmp );
#if defined(PRECISION_z) || defined(PRECISION_c)
            ptr++;
            tmp = fabs(*ptr);
            sumsq_update( 1., tmpScale, tmpSumsq, &tmp );
#endif
        }
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamRowwise of CORE_zgessq()
 */
static inline int
CORE_zgessq_row( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq )
{
    int i, j;
    double tmp;
    double *ptr, *tmpScale, *tmpSumsq;

    for(j=0; j<N; j++) {
        ptr = (double*) ( A + j * LDA );
        tmpScale = sclssq;
        tmpSumsq = sclssq+1;
        for(i=0; i<M; i++, ptr++, tmpScale+=2, tmpSumsq+=2) {
            tmp = fabs(*ptr);
            sumsq_update( 1., tmpScale, tmpSumsq, &tmp );
#if defined(PRECISION_z) || defined(PRECISION_c)
            ptr++;
            tmp = fabs(*ptr);
            sumsq_update( 1., tmpScale, tmpSumsq, &tmp );
#endif
        }
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamEltwise of CORE_zgessq()
 */
static inline int
CORE_zgessq_elt( int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq )
{
    int i, j;
    double tmp;
    double *ptr;

    for(j=0; j<N; j++) {
        ptr = (double*) ( A + j * LDA );
        for(i=0; i<M; i++, ptr++) {
            tmp = fabs(*ptr);
            sumsq_update( 1., sclssq, sclssq+1, &tmp );
#if defined(PRECISION_z) || defined(PRECISION_c)
            ptr++;
            tmp = fabs(*ptr);
            sumsq_update( 1., sclssq, sclssq+1, &tmp );
#endif
        }
    }
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgessq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to M-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*m*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
 *           i
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile A.
 * See also LAPACK zlassq.f
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = ChamColumnwise: Computes the sum of squares on each column
 *          = ChamRowwise:    Computes the sum of squares on each row
 *          = ChamEltwise:    Computes the sum of squares on all the matrix
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
 *  @param[in,out] sclssq
 *          Dimension:  (2,K)
 *          if storev == ChamColumnwise, K = N
 *          if storev == ChamRowwise,    K = M
 *          if storev == ChamEltwise,    K = 1
 *          On entry, sclssq contains K couples (sclssq[2*i], sclssq[2*i+1])
 *          which corresponds to (scale, sumsq) in the equation below
 *          ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *          respectively for the columns, the rows and the full matrix
 *          On exit, each couple is overwritten with the final result (scl, ssq).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval -k, the k-th argument had an illegal value
 *
 */
int CORE_zgessq( cham_store_t storev, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq )
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
    for (i=0; i<2*K; i+=2) {
        if ( ( sclssq[i] == -1. ) && ( sclssq[i+1] == -1. ) ) {
            sclssq[i] = 1.;
            sclssq[i+1] = 0.;
        }
    }

    if (storev == ChamColumnwise) {
        CORE_zgessq_col( M, N, A, LDA, sclssq );
    }
    else if (storev == ChamRowwise) {
        CORE_zgessq_row( M, N, A, LDA, sclssq );
    }
    else {
        CORE_zgessq_elt( M, N, A, LDA, sclssq );
    }
    return CHAMELEON_SUCCESS;
}
