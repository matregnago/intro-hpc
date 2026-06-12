/**
 *
 * @file core_zhessq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zhessq CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zhessq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to N-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*n*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( A( i, j ) ) ), abs( aimag( A( i, j ) ) ) ),
 *          i,j
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile triangular part of the
 * symmetric tile A defined by uplo.
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
 *  @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = ChamLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = ChamUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 *  @param[in] N
 *          The number of columns and rows in the tile A.
 *
 *  @param[in] A
 *          The N-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,N).
 *
 *  @param[in,out] sclssq
 *          Dimension:  (2,K)
 *          if storev == ChamColumnwise, K = N
 *          if storev == ChamRowwise,    K = N
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
int CORE_zhessq( cham_store_t storev, cham_uplo_t uplo, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *sclssq )
{
    return CORE_zsyssq( storev, uplo, N, A, LDA, sclssq );
}
