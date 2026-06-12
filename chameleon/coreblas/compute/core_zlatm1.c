/**
 *
 * @file core_zlatm1.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_dlatm1 CPU kernel
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"
#include "coreblas/random.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Double_t
 *
 * @brief Computes the entries of D(1..N) as specified by
 *    MODE, COND and IRSIGN, and as specifief in LAPACK dlatm1.f
 *
 * @warning This prototype may evolve
 *
 *******************************************************************************
 *
 * @param[in] MODE
 *           On entry describes how D is to be computed:
 *           MODE = 0 means do not change D.
 *           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
 *           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
 *           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
 *           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
 *           MODE = 5 sets D to random numbers in the range
 *                    ( 1/COND , 1 ) such that their logarithms
 *                    are uniformly distributed.
 *           MODE = 6 set D to random numbers from same distribution
 *                    as the rest of the matrix.
 *           MODE < 0 has the same meaning as ABS(MODE), except that
 *              the order of the elements of D is reversed.
 *           Thus if MODE is positive, D has entries ranging from
 *              1 to 1/COND, if negative, from 1/COND to 1,
 *           Not modified.
 *
 * @param[in] COND
 *          COND is DOUBLE PRECISION
 *           On entry, used as described under MODE above.
 *           If used, it must be >= 1. Not modified.
 *
 * @param[in] IRSIGN
 *          IRSIGN is INTEGER
 *           On entry, if MODE neither -6, 0 nor 6, determines sign of
 *           entries of D
 *           0 => leave entries of D unchanged
 *           1 => multiply each entry of D by 1 or -1 with probability .5
 *
 * @param[in] DIST
 *          DIST is INTEGER
 *           On entry, DIST specifies the type of distribution to be
 *           used to generate a random matrix. Only for MODE 6 and -6.
 *           ChamDistUniform   => UNIFORM( 0, 1 )
 *           ChamDistSymmetric => UNIFORM( -1, 1 )
 *           ChamDistNormal    => NORMAL( 0, 1 )
 *           ChamDistNormal is not supported for now.
 *
 * @param[in] ISEED
 *           Specifies the seed of the random number
 *           generator (see coreblas/random.h).
 *
 * @param[in,out] D
 *          D is DOUBLE PRECISION array, dimension ( N )
 *           Array to be computed according to MODE, COND and IRSIGN.
 *           May be changed on exit if MODE is nonzero.
 *
 * @param[in] N
 *          N is INTEGER
 *           Number of entries of D. Not modified.
 *
 *******************************************************************************
 *
 *  @retval  0  => normal termination
 *  @retval -1  => if MODE not in range -6 to 6
 *  @retval -2  => if MODE neither -6, 0 nor 6, and IRSIGN neither 0 nor 1
 *  @retval -3  => if MODE neither -6, 0 nor 6 and COND less than 1
 *  @retval -4  => if MODE equals 6 or -6 and DIST not in range 1 to 3
 *  @retval -7  => if N negative
 *
 */
int CORE_zlatm1( int MODE, double COND, int IRSIGN, cham_dist_t DIST,
                 unsigned long long int seed,
                 CHAMELEON_Complex64_t *D, int N )
{
    unsigned long long int ran;
    int i;
    double alpha, temp;
    CHAMELEON_Complex64_t ctemp;

    /*
     * Quick return if possible
     */
    if ( N == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    /*
     * Check parameters
     */
    if ( (MODE < -6) || (MODE > 6) ) {
        coreblas_error(1, "illegal value of MODE");
        return -1;
    }
    if ( ( (MODE != -6) && (MODE != 0) && (MODE != 6) ) &&
         ( (IRSIGN != 0) && (IRSIGN != 1) ) )
    {
        coreblas_error(2, "illegal value of IRSIGN");
        return -2;
    }
    if ( ( (MODE != -6) && (MODE != 0) && (MODE != 6) ) &&
         ( COND < 1. ) )
    {
        coreblas_error(3, "illegal value of COND");
        return -3;
    }
    if ( ( (MODE == -6) || (MODE == 6) ) &&
         ( (DIST < ChamDistUniform) || (DIST > ChamDistSymmetric) ) )
    {
        coreblas_error(4, "illegal value of DIST");
        return -4;
    }
    if ( N < 0 ) {
        coreblas_error(7, "illegal value of N");
        return -7;
    }

    /* Initialize the random value */
    ran = CORE_rnd64_jump( 1, seed );

    switch( abs(MODE) ) {
        /*
         * One large D value
         */
    case 1:
        D[0] = 1.;
        alpha = 1. / COND;
        for( i=1; i<N; i++ ) {
            D[i] = alpha;
        }
        break;

        /*
         * One small D value
         */
    case 2:
        for( i=0; i<N-1; i++ ) {
            D[i] = 1.;
        }
        D[N-1] = 1. / COND;
        break;

        /*
         * Exponentially distributed D values:
         */
    case 3:
        D[0] = 1.;
        if ( N > 1 ) {
            alpha = pow( COND, -1. / ((double)(N-1)) );

            for( i=1; i<N; i++ ) {
                D[i] = pow( alpha, i );
            }
        }
        break;

        /*
         * Arithmetically distributed D values
         */
    case 4:
        D[0] = 1.;
        if ( N > 1 ) {
            temp = 1. / COND;
            alpha = ( 1. - temp ) / ((double)(N-1));

            for( i=1; i<N; i++ ) {
                D[i] = ((double)(N-i-1)) * alpha + temp;
            }
        }
        break;

        /*
         * Randomly distributed D values on ( 1/COND , 1)
         */
    case 5:
        alpha = log( 1. / COND );
        for( i=0; i<N; i++ ) {
            temp = CORE_dlaran( &ran );
            D[i] = exp( alpha * temp );
        }
        break;

        /*
         * Randomly distributed D values from DIST
         */
    case 6:
        CORE_zplrnt( N, 1, D, N, N, 1, 0, seed );

        if ( DIST == ChamDistSymmetric ) {
            for( i=0; i<N; i++ ) {
                D[i] = 2. * D[i] - 1.;
            }
        }
        break;
    case 0:
    default:
        /* Nothing to be done */
        ;
    }

    /*
     * If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
     * random signs to D
     */
    if ( ( (abs(MODE) != -6) && (MODE != 0) ) &&
         (IRSIGN == 1) )
    {
        for( i=0; i<N; i++ ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
            double t1 = CORE_dlaran( &ran );
            double t2 = CORE_dlaran( &ran );
            ctemp = sqrt( -2 * log( t1 ) ) * exp( I * 2. * M_PI * t2 );
            D[i] = D[i] * ( ctemp / cabs(ctemp) );
#else
            if ( CORE_dlaran( &ran ) > .5 ) {
                D[i] = -D[i];
            }
#endif
        }
    }

    /*
     * Reverse if MODE < 0
     */
    if ( MODE < 0 ) {
        for( i=0; i<N/2; i++ ) {
            ctemp    = D[i];
            D[i]     = D[N-1-i];
            D[N-1-i] = ctemp;
        }
    }

    return CHAMELEON_SUCCESS;
}
