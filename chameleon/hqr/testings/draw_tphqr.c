/**
 *
 * @file draw_tphqr.c
 *
 * Binary to draw hierarchical trees used to compute the tsqrt/tslqt of a full matrix.
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-28
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libhqr.h>
#include "common.h"

int
main(int argc, char ** argv)
{
    libhqr_tree_t qrtree;
    libhqr_matrix_t matrix;
    int iparam[IPARAM_SIZEOF];
    int rc;

    /* Get options */
    parse_arguments( &argc, &argv, iparam );
    PASTE_CODE_IPARAM_LOCALS( iparam );

    matrix.nodes = nodes;
    matrix.p     = P;
    matrix.mt    = MT;
    matrix.nt    = NT;

    rc = libhqr_init_tphqr( &qrtree, LIBHQR_TSQR, NT, NT-1, &matrix, LIBHQR_FLAT_TREE, qr_a, qr_p );
    if ( rc != 0 ) {
        return EXIT_FAILURE;
    }
    libhqr_print_svg( &qrtree, "tphqr.svg" );
    libhqr_finalize( &qrtree );

    rc = libhqr_init_tshqr( &qrtree, 0, LIBHQR_TSQR, NT, NT-1, &matrix );
    if ( rc != 0 ) {
        return EXIT_FAILURE;
    }
    qrtree.facto = LIBHQR_QR;
    libhqr_print_svg( &qrtree, "tshqr_flat.svg" );
    libhqr_finalize( &qrtree );

    rc = libhqr_init_tshqr( &qrtree, 1, LIBHQR_TSQR, NT, NT-1, &matrix );
    if ( rc != 0 ) {
        return EXIT_FAILURE;
    }
    qrtree.facto = LIBHQR_QR;
    libhqr_print_svg( &qrtree, "tshqr_rr.svg" );
    libhqr_finalize( &qrtree );
    return EXIT_SUCCESS;
}
