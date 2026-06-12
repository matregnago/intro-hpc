/**
 *
 * @file low_fibonacci.c
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
 * Functions for low level fibonacci tree
 *
 */
#include "libhqr_internal.h"
#include <stdlib.h>

/* Return the pivot to use for the row m at step k */
static inline int
hqr_low_fibonacci_currpiv( const hqr_subpiv_t *qrpiv, int k, int m ) {
    int k_a = qrpiv->domino ? k / qrpiv->a : (k + qrpiv->p - 1 - m%(qrpiv->p)) / qrpiv->p / qrpiv->a;
    return (qrpiv->ipiv)[ k_a * (qrpiv->ldd) + ( (m / qrpiv->p) / qrpiv->a ) ];
}

/* Return the last row which has used the row m as a pivot in step k before the row start */
static inline int
hqr_low_fibonacci_prevpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
    int i;
    int k_a = qrpiv->domino ? k / qrpiv->a :  (k + qrpiv->p - 1 - p%(qrpiv->p)) / qrpiv->p / qrpiv->a;
    int p_pa = (p / qrpiv->p ) / qrpiv->a;

    for( i=start_pa+1; i<(qrpiv->ldd); i++ ) {
        if ( (qrpiv->ipiv)[i +  k_a * (qrpiv->ldd)] == p_pa ) {
            return i;
        }
    }
    return i;
}

/* Return the next row which will use the row m as a pivot in step k after it has been used by row start */
static inline int
hqr_low_fibonacci_nextpiv( const hqr_subpiv_t *qrpiv, int k, int p, int start_pa ) {
    int i;
    int k_a = qrpiv->domino ? k / qrpiv->a :  (k + qrpiv->p - 1 - p%(qrpiv->p)) / qrpiv->p / qrpiv->a;
    int p_pa = (p / qrpiv->p ) / qrpiv->a;

    for( i=start_pa-1; i>k_a; i-- ) {
        if ( (qrpiv->ipiv)[i + k_a * (qrpiv->ldd)] == p_pa ) {
            return i;
        }
    }
    return qrpiv->ldd;
}

/**
 * @brief Initialize the function pointers for a low-level fibonacci tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a fibonacci low level tree for QR/LQ factorization
 * @param[in] minMN
 *             The number of step to perform.
 */
void
hqr_low_fibonacci_init(hqr_subpiv_t *arg, int minMN) {
    int *ipiv;
    int mt;

    arg->currpiv = hqr_low_fibonacci_currpiv;
    arg->nextpiv = hqr_low_fibonacci_nextpiv;
    arg->prevpiv = hqr_low_fibonacci_prevpiv;

    mt = arg->ldd;

    arg->ipiv = (int*)calloc( mt * minMN, sizeof(int) );
    ipiv = arg->ipiv;

    /*
     * Fibonacci of order 1:
     *    u_(n+1) = u_(n) + 1
     */
    {
        int f1, k, m;

        /* Fill in the first column */
        f1 = 1;
        for (m=1; m < mt; ) {
            for (k=0; (k < f1) && (m < mt); k++, m++) {
                ipiv[m] = m - f1;
            }
            f1++;
        }

        for( k=1; k<minMN; k++) {
            for(m=k+1; m < mt; m++) {
                ipiv[ k * mt + m ] = ipiv[ (k-1) * mt + m - 1 ] + 1;
            }
        }
    }
}
