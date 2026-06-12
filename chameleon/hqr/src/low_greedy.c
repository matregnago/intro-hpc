/**
 *
 * @file low_greedy.c
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
 * Functions for low level greedy tree
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>
#include <string.h>

/* Return the pivot to use for the row m at step k */
static inline int
hqr_low_greedy_currpiv( const hqr_subpiv_t *qrpiv, int k, int m ) {
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
hqr_low_greedy_prevpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
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
hqr_low_greedy_nextpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
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
 * @brief Initialize the function pointers for a low-level greedy tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a greedy low level tree for QR/LQ factorization
 * @param[in] minMN
 *             The number of step to perform.
 */
void
hqr_low_greedy_init(hqr_subpiv_t *arg, int minMN) {
    int *ipiv;
    int mt, a, p, pa, domino;

    arg->currpiv = hqr_low_greedy_currpiv;
    arg->nextpiv = hqr_low_greedy_nextpiv;
    arg->prevpiv = hqr_low_greedy_prevpiv;

    mt = arg->ldd;
    a = arg->a;
    p = arg->p;
    pa = p * a;
    domino = arg->domino;

    if ( domino )
    {
        int j, k, start, end, firstk = 0;
        int *nT, *nZ;

        arg->minMN =  libhqr_imin( minMN, mt*a );
        minMN = arg->minMN;

        arg->ipiv = (int*)calloc( mt * minMN, sizeof(int) );
        ipiv = arg->ipiv;

        nT = (int*)calloc(minMN, sizeof(int));
        nZ = (int*)calloc(minMN, sizeof(int));

        nT[0] = mt;

        k = 0;
        while ( (!( ( nT[minMN-1] == mt - ( (minMN - 1) / a ) ) &&
                    ( nZ[minMN-1]+1 == nT[minMN-1] ) ) )
                && ( firstk < minMN ) )
        {
            int height = (nT[k] - nZ[k]) / 2;
            if ( height == 0 ) {
                while ( (firstk < minMN) &&
                        ( nT[firstk] == mt - ( firstk / a ) ) &&
                        ( nZ[firstk]+1 == nT[firstk] ) )
                {
                    if (  (( firstk % a) != a-1 )
                          && ( firstk < minMN-1 ) )
                    {
                        nT[firstk+1]++;
                    }
                    firstk++;
                }
                k = firstk;
                continue;
            }

            if (k < minMN-1) {
                nT[k+1] += height;
            }
            start = mt - nZ[k] - 1;
            end = start - height;
            nZ[k] += height;

            for( j=start; j > end; j-- ) {
                ipiv[ k*mt + j ] = (j - height);
            }

            k++;
            if (k > minMN-1) {
                k = firstk;
            }
        }

        free(nT);
        free(nZ);
    }
    else
    {
        int j, k, myrank, height, start, end;
        int *nT2DO = (int*)malloc(minMN*sizeof(int));
        int *nT = (int*)malloc(minMN*sizeof(int));
        int *nZ = (int*)malloc(minMN*sizeof(int));

        arg->ipiv = (int*)calloc( mt * minMN * p, sizeof(int) );
        ipiv = arg->ipiv;

        for ( myrank=0; myrank<p; myrank++ ) {
            int firstk;
            int lminMN = minMN;

            memset( nT2DO, 0, minMN*sizeof(int));
            memset( nT,    0, minMN*sizeof(int));
            memset( nZ,    0, minMN*sizeof(int));

            nT[0] = mt;

            for(k=0; k<lminMN; k++) {
                nT2DO[k] = libhqr_imax( mt - ((k + p - 1 - myrank) / pa), 0 );
                if ( nT2DO[k] == 0 ) {
                    lminMN = k;
                    break;
                }
            }

            k = 0;
            firstk = 0;
            while ( (!( ( nT[lminMN-1] == nT2DO[lminMN-1] ) &&
                        ( nZ[lminMN-1]+1 == nT[lminMN-1] ) ) )
                    && ( firstk < lminMN ) ) {
                height = (nT[k] - nZ[k]) / 2;
                if ( height == 0 ) {
                    while ( (firstk < lminMN) &&
                            ( nT[firstk] == nT2DO[firstk] ) &&
                            ( nZ[firstk]+1 == nT[firstk] ) )
                    {
                        if (  ( firstk < lminMN-1 )  &&
                              (( (firstk) % pa) != ((a-1)*p+myrank) ) )
                        {
                            nT[firstk+1]++;
                        }
                        firstk++;
                    }
                    k = firstk;
                    continue;
                }

                if (k < lminMN-1) {
                    nT[k+1] += height;
                }
                start = mt - nZ[k] - 1;
                end = start - height;
                nZ[k] += height;

                for( j=start; j > end; j-- ) {
                    ipiv[ myrank*mt*minMN + k*mt + j ] = (j - height);
                }

                k++;
                if (k > lminMN-1) {
                    k = firstk;
                }
            }
        }

        free(nT2DO);
        free(nT);
        free(nZ);
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
