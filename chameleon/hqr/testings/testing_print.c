/**
 *
 * @file testing_print.c
 *
 * Testing file to check the print functions.
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "libhqr.h"
#include "common.h"

int
main(int argc, char ** argv)
{
    libhqr_tree_t qrtree;
    libhqr_matrix_t matrix;

    int iparam[IPARAM_SIZEOF];
    int k, minMN, rc;

    /* Get options */
    parse_arguments( &argc, &argv, iparam );
    PASTE_CODE_IPARAM_LOCALS( iparam );

    matrix.nodes = 1;
    matrix.p = 1;
    matrix.mt = MT;
    matrix.nt = NT;

    minMN = ( NT < MT ) ? NT : MT;
    rc = libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix,
                          llvl, hlvl, qr_a, qr_p, domino, tsrr );
    if ( rc != 0 ) {
        return EXIT_FAILURE;
    }

    libhqr_print_type( &qrtree );
    libhqr_print_pivot( &qrtree );
    libhqr_print_nbgeqrt( &qrtree );

    for(k=0; k<minMN; k++) {
        libhqr_print_geqrt_k( &qrtree, k );
        libhqr_print_next_k( &qrtree, k );
        libhqr_print_prev_k( &qrtree, k );
    }

    libhqr_finalize( &qrtree );

    return EXIT_SUCCESS;
}
