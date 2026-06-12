/**
 *
 * @file low_adaptiv.c
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 * Functions for low level adaptiv tree for SVD/EVD reduction to band.
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

/* Return the pivot to use for the row m at step k */
inline static int
svd_low_adaptiv_currpiv( const hqr_subpiv_t *arg,
                         int k, int m )
{
    int *ipiv = arg->ipiv + (m%arg->p * arg->minMN + k) * arg->ldd;
    int a = ipiv[0];

    ipiv+=2;
    return ipiv[ ( m  / arg->p ) / a ];
}

/* Return the last row which has used the row m as a pivot in step k before the row start */
inline static int
svd_low_adaptiv_prevpiv( const hqr_subpiv_t *arg,
                         int k, int p, int start_pa )
{
    int i;
    int *ipiv = arg->ipiv + (p%arg->p * arg->minMN + k) * arg->ldd;
    int a = ipiv[0];
    int ldd = ipiv[1];
    int p_pa = p / arg->p / a;

    ipiv+=2;
    for( i=start_pa+1; i<ldd; i++ ) {
        if ( ipiv[i] == p_pa ) {
            return i;
        }
    }
    return i;
}

/* Return the next row which will use the row m as a pivot in step k after it has been used by row start */
inline static int
svd_low_adaptiv_nextpiv( const hqr_subpiv_t *arg,
                         int k, int p, int start_pa )
{
    int i;
    int *ipiv = arg->ipiv + (p%arg->p * arg->minMN + k ) * arg->ldd;
    int a = ipiv[0];
    int ldd = ipiv[1];
    int pa = arg->p * a;
    int k_a = (k + arg->p - 1 - p%(arg->p)) / arg->p / a;
    int p_pa = p / pa;

    ipiv+=2;
    for( i=start_pa-1; i> k_a; i-- ) {
        if ( ipiv[i] == p_pa ) {
            return i;
        }
    }

    return ldd;
}

/**
 * @brief Initialize the function pointers for a low-level adaptiv tree for SVD computations
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define an adaptive low level tree in the general to band
 *             reduction process of the SVD computation.
 * @param[in] gmt
 *             The global number of rows of tiles in the matrix to consider.
 * @param[in] gnt
 *             The global number of columns of tiles in the matrix to consider.
 * @param[in] nbcores
 *             The number of cores available on each node to perfrom the computation.
 * @param[in] ratio
 *             The number of tiles per thread that should be available at each
 *             step of the factorization.
 */
void
svd_low_adaptiv_init( hqr_subpiv_t *arg,
                      int gmt, int gnt, int nbcores, int ratio )
{
    int *ipiv;
    int mt, a, p, pa, maxmt, myrank;
    int j, k, height, start, end, nT, nZ;
    int minMN = libhqr_imin(gmt, gnt);

    arg->currpiv = svd_low_adaptiv_currpiv;
    arg->nextpiv = svd_low_adaptiv_nextpiv;
    arg->prevpiv = svd_low_adaptiv_prevpiv;

    p = arg->p;

    end = 0;

    /**
     * Compute the local greedy tree of each column, on each node
     */
    maxmt = 1;
    for(k=0; k<minMN; k++) {
        /**
         * The objective is to have at least two columns of TS to reduce per
         * core, so it must answer the following inequality:
         * ((gmt-k) / p / a ) * (gnt-k) >= ( ratio * nbcores );
         * so,
         * a <= mt * (gnt-k) / (ratio * nbcores )
         */
        height = libhqr_iceil( gmt-k, p );
        a = libhqr_imax( height * (gnt-k) / (ratio * nbcores), 1 );

        /* Now let's make sure all sub-parts are equilibrate */
        j = libhqr_iceil( height, a );
        a = libhqr_iceil( gmt-k, j );

        /* Compute max dimension of the tree */
        mt = libhqr_iceil( gmt, p * a );
        maxmt = libhqr_imax( mt, maxmt );
    }

    arg->ldd = maxmt + 2;
    arg->ipiv = (int*)calloc( arg->ldd * minMN * p, sizeof(int) );
    ipiv = arg->ipiv;

    for ( myrank=0; myrank<p; myrank++ ) {

        /**
         * Compute the local greedy tree of each column, on each node
         */
        for(k=0; k<minMN; k++, ipiv += arg->ldd) {
            /**
             * The objective is to have at least two columns of TS to reduce per
             * core, so it must answer the following inequality:
             * (ldd / a ) * (gnt-k) >= ( ratio * nbcores );
             * so,
             * a <= mt * (gnt-k) / (ratio * nbcores )
             */
            height = libhqr_iceil( gmt-k, p );
            a = libhqr_imax( height * (gnt-k) / (ratio * nbcores), 1 );

            /* Now let's make sure all sub-parts are equilibrate */
            j = libhqr_iceil( height, a );
            a = libhqr_iceil( gmt-k, j );

            pa = p * a;
            mt = libhqr_iceil( gmt, pa );
            ipiv[0] = a;
            ipiv[1] = mt;

            assert( a  > 0 );
            assert( mt < arg->ldd-1 );

            /* Number of tiles to factorized in this column on this rank */
            nT = libhqr_imax( mt - ((k + p - 1 - myrank) / pa), 0 );
            /* Number of tiles already killed */
            nZ = 0;

            assert( nT <= mt );

            /* No more computations on this node */
            if ( nT == 0 ) {
                continue;
            }

            while( nZ < (nT-1) ) {
                height = (nT - nZ) / 2;
                start = mt - nZ - 1;
                end = start - height;
                nZ += height;

                for( j=start; j > end; j-- ) {
                    ipiv[ j+2 ] = (j - height);
                }
            }
            assert( nZ+1 == nT );
        }
    }

#if 0
    {
        int m, k;
        for(m=0; m<mt; m++) {
            printf("%3d | ", m);
            for (k=0; k<minMN; k++) {
                printf( "%3d ", ipiv[k*(arg->ldd) + m] );
            }
            printf("\n");
        }
    }
    if (!arg->domino) {
        int m, k, myrank;
        for ( myrank=1; myrank<p; myrank++ ) {
            ipiv += arg->ldd * minMN;
            printf("-------- rank %d ---------\n", myrank );
            for(m=0; m<mt; m++) {
                printf("%3d | ", m);
                for (k=0; k<minMN; k++) {
                    int k_a = (k + p - 1 - myrank) / p / a;
                    if ( m >= k_a ) {
                        printf( "%3d ", ipiv[k * arg->ldd + m] );
                    }
                    else {
                        printf( "--- " );
                    }
                }
                printf("\n");
            }
        }
    }
#endif
}
