/**
 *
 * @file high_fibonacci.c
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
 * Functions for high level fibonacci tree, and init for duplicated greedy.
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

/****************************************************
 *          HQR_HIGH_FIBONACCI_TREE
 ***************************************************/
/* Return the pivot to use for the row m at step k */
static inline int
hqr_high_fibonacci_currpiv( const hqr_subpiv_t *qrpiv, int k, int m ) {
    return (qrpiv->ipiv)[ m-k ] + k;
}

/* Return the last row which has used the row m as a pivot in step k before the row start */
static inline int
hqr_high_fibonacci_prevpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start ) {
    int i;
    myassert( p >= k && start >= p && start-k <= qrpiv->p);

    int lp    = p - k;
    int lstart= start - k;
    int end   = libhqr_imin(qrpiv->ldd-k, qrpiv->p);
    for( i=lstart+1; i<end; i++ ) {
        if ( (qrpiv->ipiv)[i] == lp ) {
            return i+k;
        }
    }
    return qrpiv->ldd;
}

/* Return the next row which will use the row m as a pivot in step k after it has been used by row start */
static inline int
hqr_high_fibonacci_nextpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start ) {
    int i;
    myassert( p>=k && (start == qrpiv->ldd || start-k <= qrpiv->p) );

    for( i=libhqr_imin(start-k-1, qrpiv->p-1); i>0; i-- ) {
        if ( (qrpiv->ipiv)[i] == (p-k) ) {
            return i + k;
        }
    }
    return qrpiv->ldd;
}

/**
 * @brief Initialize the function pointers for a high-level fibonacci tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a fibonacci high level tree for QR/LQ factorization.
 */
void
hqr_high_fibonacci_init(hqr_subpiv_t *arg) {
    int *ipiv;
    int p;

    arg->currpiv = hqr_high_fibonacci_currpiv;
    arg->nextpiv = hqr_high_fibonacci_nextpiv;
    arg->prevpiv = hqr_high_fibonacci_prevpiv;

    p = arg->p;

    arg->ipiv = (int*)calloc( p, sizeof(int) );
    ipiv = arg->ipiv;

    /*
     * Fibonacci of order 1:
     *    u_(n+1) = u_(n) + 1
     */
    {
        int f1, k, m;

        /* Fill in the first column */
        f1 = 1;
        for (m=1; m < p; ) {
            for (k=0; (k < f1) && (m < p); k++, m++) {
                ipiv[m] = m - f1;
            }
            f1++;
        }
    }
}

/**
 * @brief Initialize the function pointers for a replicated high-level greedy tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a replicated greedy high level tree for QR/LQ factorization.
 *             The greedy distribution is perfoned once, and then replicated on
 *             all other columns.
 */
void hqr_high_greedy1p_init(hqr_subpiv_t *arg){
    int *ipiv;
    int mt, p;

    arg->currpiv = hqr_high_fibonacci_currpiv;
    arg->nextpiv = hqr_high_fibonacci_nextpiv;
    arg->prevpiv = hqr_high_fibonacci_prevpiv;

    mt = arg->ldd;
    p = arg->p;

    arg->ipiv = (int*)calloc( p, sizeof(int) );
    ipiv = arg->ipiv;

    {
        int j, height, start, end;
        int nT, nZ;

        nT = mt;
        nZ = libhqr_imax( mt - p, 0 );
        while ( !( ( nT == mt ) &&
                   ( nT == nZ+1 ) ) )
        {
            height = (nT - nZ) / 2;
            if ( height == 0 ) {
                break;
            }

            start = mt - nZ - 1;
            end = start - height;
            nZ += height;

            for( j=start; j > end; j-- ) {
                ipiv[ j ] = (j - height);
            }
        }
    }
}
