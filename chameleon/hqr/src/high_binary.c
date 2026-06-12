/**
 *
 * @file high_binary.c
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
 * Functions for high level binary tree
 *
 */
#include "libhqr_internal.h"
#include <math.h>

static int
hqr_high_binary_currpiv(const hqr_subpiv_t *arg, int k, int m)
{
    int tmp1 = m - k;
    int tmp2 = 1;
    (void)arg;

    if ( tmp1 == 0) {
        return 0;
    }
    while( (tmp1 != 0 ) && (tmp1 % 2 == 0) ) {
        tmp1 = tmp1 >> 1;
        tmp2 = tmp2 << 1;
    }
    assert( m - tmp2 >= k );
    return m - tmp2;
}

static int
hqr_high_binary_nextpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    int tmpp, bit;
    myassert( (start == arg->ldd) || (hqr_high_binary_currpiv( arg, k, start ) == p) );

    if ( start <= p ) {
        return arg->ldd;
    }

    int offset = p - k;
    bit = 0;
    if (start != arg->ldd) {
        while( ( (start - k) & (1 << bit ) ) == 0 ) {
            bit++;
        }
        bit++;
    }

    tmpp = offset | (1 << bit);
    if ( ( tmpp != offset ) && (tmpp < arg->p) && ( tmpp+k < arg->ldd ) ) {
        return tmpp + k;
    }
    else {
        return arg->ldd;
    }
}

static int
hqr_high_binary_prevpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    int offset = p - k;

    myassert( start >= p && ( start == p || hqr_high_binary_currpiv( arg, k, start ) == p ) );

    if ( (start == p) && ( offset%2 == 0 ) ) {
        int i, bit;
        if ( offset == 0 ) {
            bit = (int)( log( (double)( libhqr_imin(arg->p, arg->ldd - k) ) ) / log( 2. ) );
        }
        else {
            bit = 0;
            while( (offset & (1 << bit )) == 0 ) {
                bit++;
            }
        }
        for( i=bit; i>-1; i--){
            int tmp = offset | (1 << i);
            if ( ( offset != tmp ) && ( tmp < arg->p ) && (tmp+k < arg->ldd) ) {
                return tmp+k;
            }
        }
        return arg->ldd;
    }

    if ( (start - p) > 1 ) {
        return p + ( (start - p) >> 1 );
    }
    else {
        return arg->ldd;
    }
}

/**
 * @brief Initialize the function pointers for a high-level binary tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a binary high level tree for QR/LQ factorization.
 */
void
hqr_high_binary_init(hqr_subpiv_t *arg) {
    arg->currpiv = hqr_high_binary_currpiv;
    arg->nextpiv = hqr_high_binary_nextpiv;
    arg->prevpiv = hqr_high_binary_prevpiv;
    arg->ipiv = NULL;
}

