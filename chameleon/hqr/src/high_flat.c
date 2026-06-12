/**
 *
 * @file high_flat.c
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
 * Functions for high level flat tree
 *
 */
#include "libhqr_internal.h"

static int
hqr_high_flat_currpiv(const hqr_subpiv_t *arg, int k, int m)
{
    (void)arg;
    (void)m;
    return k;
}

static int
hqr_high_flat_nextpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    if ( p == k && arg->ldd > 1 ) {
        if ( start == arg->ldd ) {
            return p+1;
        }
        else if ( start < arg->ldd && (start-k < arg->p-1) ) {
            return start+1;
        }
    }
    return arg->ldd;
}

static int
hqr_high_flat_prevpiv(const hqr_subpiv_t *arg, int k, int p, int start)
{
    assert( arg->p > 1 );
    if ( p == k && arg->ldd > 1 ) {
        if ( start == p && p != arg->ldd-1 ) {
            return libhqr_imin( p + arg->p - 1, arg->ldd - 1 );
        }
        else if ( start > p + 1 && (start-k < arg->p)) {
            return start-1;
        }
    }
    return arg->ldd;
}

/**
 * @brief Initialize the function pointers for a high-level flat tree
 *
 * @param[in,out] arg
 *             The data structure to initialize. On exit, it can be used to
 *             define a flat high level tree for QR/LQ factorization.
 */
void
hqr_high_flat_init(hqr_subpiv_t *arg) {
    arg->currpiv = hqr_high_flat_currpiv;
    arg->nextpiv = hqr_high_flat_nextpiv;
    arg->prevpiv = hqr_high_flat_prevpiv;
    arg->ipiv = NULL;
}

