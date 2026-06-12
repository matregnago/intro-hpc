/**
 *
 * @file high_greedy.c
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
 * Functions for high level greedy tree
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

static int
hqr_high_greedy_currpiv(const hqr_subpiv_t *arg, int k, int m)
{
    myassert( m >= k && m < k+arg->p );
    return (arg->ipiv)[ k * (arg->p) + (m - k) ];
}

static int
hqr_high_greedy_nextpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    int i;
    myassert( (start >= k && start < k+arg->p) || start == arg->ldd );
    for( i=libhqr_imin(start-1, k+arg->p-1); i > k; i-- ) {
        if ( (arg->ipiv)[i-k + k* (arg->p)] == p ) {
            return i;
        }
    }
    return arg->ldd;
}

static int
hqr_high_greedy_prevpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    int i;
    myassert( (start >= k && start < k+arg->p) || start == p );
    for( i=start-k+1; i<arg->p; i++ ) {
        if ( (arg->ipiv)[i +  k * (arg->p)] == p ) {
            return k+i;
        }
    }
    return arg->ldd;
}

/**
 * @brief Initialize the function pointers for a high-level greedy tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a greedy high level tree for QR/LQ factorization.
 *             The greedy algorithm is applied on all columns taking into
 *             account the history of previous steps.
 */
void
hqr_high_greedy_init(hqr_subpiv_t *arg, int minMN) {
    int *ipiv;
    int mt, p;

    arg->currpiv = hqr_high_greedy_currpiv;
    arg->nextpiv = hqr_high_greedy_nextpiv;
    arg->prevpiv = hqr_high_greedy_prevpiv;

    mt = arg->ldd;
    p = arg->p;

    arg->ipiv = (int*)calloc( p * minMN, sizeof(int) );
    ipiv = arg->ipiv;

    {
        int j, k, height, start, end, firstk = 0;
        int *nT = (int*)calloc(minMN, sizeof(int));
        int *nZ = (int*)calloc(minMN, sizeof(int));

        nT[0] = mt;
        nZ[0] = libhqr_imax( mt - p, 0 );
        for(k=1; k<minMN; k++) {
            height = libhqr_imax(mt-k-p, 0);
            nT[k] = height;
            nZ[k] = height;
        }

        k = 0;
        while ( (!( ( nT[minMN-1] == mt - (minMN - 1) ) &&
                    ( nZ[minMN-1]+1 == nT[minMN-1] ) ) )
                && ( firstk < minMN ) ) {
            height = (nT[k] - nZ[k]) / 2;
            if ( height == 0 ) {
                while ( (firstk < minMN) &&
                        ( nT[firstk] == mt - firstk ) &&
                        ( nZ[firstk]+1 == nT[firstk] ) ) {
                    firstk++;
                }
                k = firstk;
                continue;
            }

            start = mt - nZ[k] - 1;
            end = start - height;
            nZ[k] += height;
            if (k < minMN-1) {
                nT[k+1] = nZ[k];
            }

            for( j=start; j > end; j-- ) {
                ipiv[ k*p + j-k ] = (j - height);
            }

            k++;
            if (k > minMN-1) {
                k = firstk;
            }
        }

        free(nT);
        free(nZ);
    }
}
