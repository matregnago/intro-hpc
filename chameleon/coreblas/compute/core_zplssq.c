/**
 *
 * @file core_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zplssq CPU kernel
 *
 * @version 1.4.0
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
 * @brief Subcase storev == ChamColumnwise of CORE_zplssq()
 */
static inline int
CORE_zplssq_col( int M, int N,
                 double *sclssqin,
                 double *sclssqout )
{
    int i, j;
    double *tmpScaleIn, *tmpSumsqIn, *tmpScaleOut, *tmpSumsqOut;

    tmpScaleIn = sclssqin;
    tmpSumsqIn = sclssqin+1;
    tmpScaleOut = sclssqout;
    tmpSumsqOut = sclssqout+1;

    for(j=0; j<N; j++) {
        for(i=0; i<M; i++) {
            sumsq_update_2( tmpScaleIn, tmpSumsqIn, tmpScaleOut, tmpSumsqOut );
            tmpScaleIn+=2;
            tmpSumsqIn+=2;
        }
        tmpScaleOut+=2;
        tmpSumsqOut+=2;
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamRowwise of CORE_zplssq()
 */
static inline int
CORE_zplssq_row( int M, int N,
                 double *sclssqin,
                 double *sclssqout )
{
    int i, j;
    double *tmpScaleIn, *tmpSumsqIn, *tmpScaleOut, *tmpSumsqOut;

    tmpScaleIn = sclssqin;
    tmpSumsqIn = sclssqin+1;

    for(j=0; j<N; j++) {
        tmpScaleOut = sclssqout;
        tmpSumsqOut = sclssqout+1;
        for(i=0; i<M; i++) {
            sumsq_update_2( tmpScaleIn, tmpSumsqIn, tmpScaleOut, tmpSumsqOut );
            tmpScaleIn+=2;
            tmpSumsqIn+=2;
            tmpScaleOut+=2;
            tmpSumsqOut+=2;
        }
    }
    return CHAMELEON_SUCCESS;
}
/**
 * @brief Subcase storev == ChamEltwise of CORE_zplssq()
 */
static inline int
CORE_zplssq_elt( int M, int N,
                 double *sclssqin,
                 double *sclssqout )
{
    int i, j;
    double *tmpScaleIn, *tmpSumsqIn, *tmpScaleOut, *tmpSumsqOut;

    tmpScaleIn = sclssqin;
    tmpSumsqIn = sclssqin+1;
    tmpScaleOut = sclssqout;
    tmpSumsqOut = sclssqout+1;

    for(j=0; j<N; j++) {
        for(i=0; i<M; i++) {
            sumsq_update_2( tmpScaleIn, tmpSumsqIn, tmpScaleOut, tmpSumsqOut );
            tmpScaleIn+=2;
            tmpSumsqIn+=2;
        }
    }
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zplssq adds the values ssq ensuring scl >= 0
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
 *          The number of rows in the tile sclssqin.
 *
 *  @param[in] N
 *          The number of columns in the tile sclssqin.
 *
 *  @param[in] sclssqin
 *          The 2*M-by-2*N matrix on which to compute the sums.
 *
 *  @param[in,out] sclssqout
 *          Dimension:  (2,K)
 *          if storev == ChamColumnwise, K = N
 *          if storev == ChamRowwise,    K = M
 *          if storev == ChamEltwise,    K = 1
 *          On entry, sclssqout contains M-by-N couples (sclssqout[2*i], sclssqout[2*i+1])
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
int CORE_zplssq( cham_store_t storev, int M, int N,
                 double *sclssqin,
                 double *sclssqout )
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
        if ( ( sclssqout[i] == -1. ) && ( sclssqout[i+1] == -1. ) ) {
            sclssqout[i] = 1.;
            sclssqout[i+1] = 0.;
        }
    }

    if (storev == ChamColumnwise) {
        CORE_zplssq_col( M, N, sclssqin, sclssqout );
    }
    else if (storev == ChamRowwise) {
        CORE_zplssq_row( M, N, sclssqin, sclssqout );
    }
    else {
        CORE_zplssq_elt( M, N, sclssqin, sclssqout );
    }
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zplssq2 computes scl*sqrt(ssq) for each couple in a vector
 *
 *******************************************************************************
 *
 *  @param[in] N
 *          The number of columns in the tile sclssq.
 *
 *  @param[in,out] sclssq
 *          The 2*N matrix on which to compute scl*sqrt(ssq)
 *          On entry contains all couple (scl,ssq) in (sclssq[i],sclssq[i+1])
 *          On exist returns scl*sqrt(ssq) stored in sclssq[2*i], i = 0, ..., N-1
 *          so that the result is stored in the first line.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval -k, the k-th argument had an illegal value
 *
 */
int CORE_zplssq2( int N,
                  double *sclssq )
{
    int i;
    for (i=0; i<2*N; i+=2) {
        sclssq[i] = sclssq[i]*sqrt(sclssq[i+1]);
    }
    return CHAMELEON_SUCCESS;
}
