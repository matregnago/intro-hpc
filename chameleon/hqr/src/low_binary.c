/**
 *
 * @file low_binary.c
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
 * Functions for low level binary tree
 *
 */
#include "libhqr_internal.h"
#include <math.h>

static inline int
hqr_low_binary_currpiv(const hqr_subpiv_t *arg, int k, int m)
{
    int k_a = arg->domino ? k / arg->a :  (k + arg->p - 1 - m%(arg->p)) / arg->p / arg->a;
    int m_pa = (m / arg->p ) / arg->a;

    int tmp1 = m_pa - k_a;
    int tmp2 = 1;
    (void)arg;

    if ( tmp1 == 0) {
        return 0;
    }
    while( (tmp1 != 0 ) && (tmp1 % 2 == 0) ) {
        tmp1 = tmp1 >> 1;
        tmp2 = tmp2 << 1;
    }
    assert( m_pa - tmp2 >= k_a );
    return m_pa - tmp2;
}

static inline int
hqr_low_binary_nextpiv(const hqr_subpiv_t *arg, int k, int p, int start_pa)
{
    int k_a = arg->domino ? k / arg->a :  (k + arg->p - 1 - p%(arg->p)) / arg->p / arg->a;
    int p_pa = (p / arg->p ) / arg->a;

    int tmpp, bit;
    myassert( (start_pa == arg->ldd) || (hqr_low_binary_currpiv( arg, k, start_pa*arg->a*arg->p ) == p_pa || !arg->domino) );

    if ( start_pa <= p_pa ) {
        return arg->ldd;
    }

    int offset = p_pa-k_a;
    bit = 0;
    if (start_pa != arg->ldd) {
        while( ( (start_pa-k_a) & (1 << bit ) ) == 0 ) {
            bit++;
        }
        bit++;
    }

    tmpp = offset | (1 << bit);
    if ( ( tmpp != offset ) && ( tmpp+k_a < arg->ldd ) ) {
        return tmpp + k_a;
    }
    else {
        return arg->ldd;
    }
}

static inline int
hqr_low_binary_prevpiv(const hqr_subpiv_t *arg, int k, int p, int start_pa)
{
    int k_a = arg->domino ? k / arg->a :  (k + arg->p - 1 - p%(arg->p)) / arg->p / arg->a;
    int p_pa = (p / arg->p ) / arg->a;
    int offset = p_pa - k_a;

    myassert( start_pa >= p_pa && ( start_pa == p_pa || !arg->domino ||
                                    hqr_low_binary_currpiv( arg, k, start_pa*arg->a*arg->p ) == p_pa ) );

    if ( (start_pa == p_pa) && ( offset%2 == 0 ) ) {
        int i, bit;
        if ((p_pa - k_a) == 0) {
            bit = (int)( log( (double)(arg->ldd - k_a) ) / log( 2. ) );
        }
        else {
            bit = 0;
            while( (offset & (1 << bit )) == 0 ) {
                bit++;
            }
        }
        for( i=bit; i>-1; i--){
            int tmp = offset | (1 << i);
            if ( ( offset != tmp ) && ( tmp+k_a < arg->ldd ) ) {
                return tmp+k_a;
            }
        }
        return arg->ldd;
    }

    if ( (start_pa - p_pa) > 1 ) {
        return p_pa + ( (start_pa - p_pa) >> 1 );
    }
    else {
        return arg->ldd;
    }
}

/**
 * @brief Initialize the function pointers for a low-level binary tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a low level binary tree.
 */
void
hqr_low_binary_init(hqr_subpiv_t *arg) {
    arg->currpiv = hqr_low_binary_currpiv;
    arg->nextpiv = hqr_low_binary_nextpiv;
    arg->prevpiv = hqr_low_binary_prevpiv;
    arg->ipiv = NULL;
}
