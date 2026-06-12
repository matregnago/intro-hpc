/**
 *
 * @file testing_systolic.c
 *
 * Testing file for all combinations of systolic trees.
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
    int rc, ret = 0;

    /* Get options */
    parse_arguments( &argc, &argv, iparam );
    PASTE_CODE_IPARAM_LOCALS( iparam );

    /* Test many combinations */
    if ( check ) {
        int all_qr_p[] = { 1, 3, 5, 7, 8 };
        int all_qr_q[] = { 1, 2, 4, 6, 7 };
        int all_mt[]   = { 1, 3, 4, 10, 17, 25, 128 };
        int all_nt[]   = { 1, 2, 5, 13, 26, 58 };

        int nb_qr_p = sizeof( all_qr_p ) / sizeof( int );
        int nb_qr_q = sizeof( all_qr_q ) / sizeof( int );
        int nb_mt   = sizeof( all_mt   ) / sizeof( int );
        int nb_nt   = sizeof( all_nt   ) / sizeof( int );

        int qrp, qrq, m, n;
        int done, todo;
        todo = 0;
        done = 0;

        todo +=  nb_mt * nb_nt * nb_qr_q * nb_qr_p;

        /* P */
        for( qrp=0; qrp<nb_qr_p; qrp++) {
            P = all_qr_p[qrp];
            matrix.p = P;

            /* Q */
            for( qrq=0; qrq<nb_qr_q; qrq++) {
                Q = all_qr_q[qrq];
                matrix.nodes = P * Q;

                /* MT */
                for( m=0; m<nb_mt; m++) {
                    MT = all_mt[m];
                    matrix.mt = MT;

                    /* NT */
                    for( n=0; n<nb_nt; n++) {
                        NT = all_nt[n];
                        matrix.nt = NT;

                        rc = libhqr_init_sys( &qrtree, LIBHQR_QR, &matrix, P, Q );
                        if ( rc == 0 ) {
                            rc = libhqr_check( &qrtree );
                            libhqr_finalize( &qrtree );
                        }

                        if (rc != 0) {
                            fprintf(stderr, "%s -M %d -N %d --qr_p=%d --qr_a=%d      FAILED(%d)\n",
                                    argv[0], MT, NT, P, Q, rc);
                            ret++;
                            if (check == 1) {
                                return 0;
                            }
                        }

                        done++;
                        printf("\r%6d / %6d", done, todo);
                    }
                }
            }
        }
    }
    else {
        matrix.nodes = nodes;
        matrix.p     = P;
        matrix.mt    = MT;
        matrix.nt    = NT;

        rc = libhqr_init_sys( &qrtree, LIBHQR_QR, &matrix,
                              qr_p, qr_a );

        if ( rc == 0 ) {
            rc = libhqr_check( &qrtree );
            libhqr_finalize( &qrtree );
        }

        if (rc != 0) {
            fprintf(stderr, "%s -M %d -N %d --qr_p=%d --qr_a=%d      FAILED(%d)\n",
                    argv[0], MT, NT, qr_p, qr_a, rc);
            ret++;
        }
        else {
            fprintf(stderr, "%s -M %d -N %d --qr_p=%d --qr_a=%d      SUCCESS\n",
                    argv[0], MT, NT, qr_p, qr_a);
        }
    }

    if ( check > 1 ) {
        printf( "<DartMeasurement name=\"failures\" type=\"numeric/integer\"\n"
                "                 encoding=\"none\" compression=\"none\">\n"
                "%d\n"
                "</DartMeasurement>\n",
                ret );
    }
    if ( ret == 0 ) {
        return EXIT_SUCCESS;
    }
    else {
        printf( "%d tests failed !!!\n", ret  );
        return EXIT_FAILURE;
    }
}
