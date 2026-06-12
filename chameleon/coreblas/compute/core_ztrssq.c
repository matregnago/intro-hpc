/**
 *
 * @file core_ztrssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztrssq CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <math.h>
#include "coreblas/lapacke.h"
#include "coreblas.h"
#include "coreblas/sumsq_update.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_ztrssq returns the values scl and ssq such that
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
 *  @param[in,out] scale
 *          On entry, the value  scale  in the equation above.
 *          On exit, scale is overwritten with the value scl.
 *
 *  @param[in,out] sumsq
 *          On entry, the value  sumsq  in the equation above.
 *          On exit, SUMSQ is overwritten with the value ssq.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval -k, the k-th argument had an illegal value
 *
 */
int CORE_ztrssq( cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *scale, double *sumsq )
{
    int i, j, imax;
    int idiag = (diag == ChamUnit) ? 1 : 0;
    double *ptr;

    if ( diag == ChamUnit ){
        double tmp = sqrt( chameleon_min(M, N) );
        sumsq_update( 1, scale, sumsq, &tmp );
    }

    if  (uplo == ChamUpper ) {
        M = chameleon_min(M, N);

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA );
            imax = chameleon_min(j+1-idiag, M);

            for(i=0; i<imax; i++, ptr++) {
		sumsq_update( 1, scale, sumsq, ptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                ptr++;
		sumsq_update( 1, scale, sumsq, ptr );
#endif
            }
        }
    }
    else {
        N = chameleon_min(M, N);

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * (LDA+1) + idiag );

            for(i=j+idiag; i<M; i++, ptr++) {
                sumsq_update( 1., scale, sumsq, ptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                ptr++;
                sumsq_update( 1., scale, sumsq, ptr );
#endif
            }
        }
    }
    return CHAMELEON_SUCCESS;
}
