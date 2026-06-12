/**
 *
 * @file low_greedy1p.c
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
 * Functions for low level duplicated greedy tree
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

/* Return the pivot to use for the row m at step k */
static inline int
hqr_low_greedy1p_currpiv( const hqr_subpiv_t *qrpiv, int k, int m ) {
    if (qrpiv->domino) {
        return (qrpiv->ipiv)[ k * (qrpiv->ldd) + ( (m / qrpiv->p) / qrpiv->a ) ];
    }
    else {
        return (qrpiv->ipiv)[ ( (m%qrpiv->p) * qrpiv->minMN + k ) * (qrpiv->ldd)
                              + ( ( m  / qrpiv->p ) / qrpiv->a ) ];
    }
}

/* Return the last row which has used the row m as a pivot in step k before the row start */
static inline int
hqr_low_greedy1p_prevpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
    int i;
    int p_pa = p / qrpiv->p / qrpiv->a;
    int *ipiv = qrpiv->domino ? qrpiv->ipiv : qrpiv->ipiv + p%qrpiv->p * qrpiv->minMN *qrpiv->ldd;

    for( i=start_pa+1; i<(qrpiv->ldd); i++ ) {
        if ( ipiv[i +  k * (qrpiv->ldd)] == p_pa ) {
            return i;
        }
    }
    return i;
}

/* Return the next row which will use the row m as a pivot in step k after it has been used by row start */
static inline int
hqr_low_greedy1p_nextpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
    int i;
    int pa = qrpiv->p * qrpiv->a;
    int k_a = qrpiv->domino ? k / qrpiv->a :  (k + qrpiv->p - 1 - p%(qrpiv->p)) / qrpiv->p / qrpiv->a;
    int p_pa = p / pa;
    int *ipiv = qrpiv->domino ? qrpiv->ipiv : qrpiv->ipiv + p%qrpiv->p * qrpiv->minMN *qrpiv->ldd;

    for( i=start_pa-1; i> k_a; i-- ) {
        if ( ipiv[i + k * (qrpiv->ldd)] == p_pa ) {
            return i;
        }
    }

    return qrpiv->ldd;
}

/**
 * @brief Initialize the function pointers for a replicated low-level greedy tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a replicated greedy low level tree for QR/LQ factorization.
 *             The greedy distribution is computed once, and replicated on all
 *             other columns.
 * @param[in] minMN
 *             The number of step to perform.
 */
void
hqr_low_greedy1p_init(hqr_subpiv_t *arg, int minMN)
{
    int *ipiv;
    int mt, a, p, pa, domino;
    int j, k, height, start, end, nT, nZ;

    arg->currpiv = hqr_low_greedy1p_currpiv;
    arg->nextpiv = hqr_low_greedy1p_nextpiv;
    arg->prevpiv = hqr_low_greedy1p_prevpiv;

    mt = arg->ldd;
    a = arg->a;
    p = arg->p;
    pa = p * a;
    domino = arg->domino;

    /* This section has not been coded yet, and will perform a classic greedy */
    if ( domino )
    {
        arg->minMN =  libhqr_imin( minMN, mt*a );
        minMN = arg->minMN;

        arg->ipiv = (int*)calloc( mt * minMN, sizeof(int) );
        ipiv = arg->ipiv;

        /**
         * Compute the local greedy tree of each column, on each node
         */
        for(k=0; k<minMN; k++) {
            /* Number of tiles to factorized in this column on this rank */
            nT = libhqr_imax( mt - (k / a), 0 );
            /* Number of tiles already killed */
            nZ = 0;

            while( nZ < (nT-1) ) {
                height = (nT - nZ) / 2;
                start = mt - nZ - 1;
                end = start - height;
                nZ += height;

                for( j=start; j > end; j-- ) {
                    ipiv[ k*mt + j ] = (j - height);
                }
            }
            assert( nZ+1 == nT );
        }
    }
    else
    {
        int myrank;
        end = 0;

        arg->ipiv = (int*)calloc( mt * minMN * p, sizeof(int) );
        ipiv = arg->ipiv;

        for ( myrank=0; myrank<p; myrank++ ) {

            /**
             * Compute the local greedy tree of each column, on each node
             */
            for(k=0; k<minMN; k++) {
                /* Number of tiles to factorized in this column on this rank */
                nT = libhqr_imax( mt - ((k + p - 1 - myrank) / pa), 0 );
                /* Number of tiles already killed */
                nZ = 0;

                /* No more computations on this node */
                if ( nT == 0 ) {
                    break;
                }

                while( nZ < (nT-1) ) {
                    height = (nT - nZ) / 2;
                    start = mt - nZ - 1;
                    end = start - height;
                    nZ += height;

                    for( j=start; j > end; j-- ) {
                        ipiv[ myrank*mt*minMN + k*mt + j ] = (j - height);
                    }
                }
                assert( nZ+1 == nT );
            }
        }
    }

#if 0
    {
        int m, k;
        for(m=0; m<mt; m++) {
            printf("%3d | ", m);
            for (k=0; k<minMN; k++) {
                printf( "%3d ", ipiv[k*mt + m] );
            }
            printf("\n");
        }
    }
    if (!arg->domino) {
        int m, k, myrank;
        for ( myrank=1; myrank<p; myrank++ ) {
            ipiv += mt*minMN;
            printf("-------- rank %d ---------\n", myrank );
            for(m=0; m<mt; m++) {
                printf("%3d | ", m);
                for (k=0; k<minMN; k++) {
                    int k_a = (k + p - 1 - myrank) / p / a;
                    if ( m >= k_a ) {
                        printf( "%3d ", ipiv[k*mt + m] );
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
