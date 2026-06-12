/**
 *
 * @file testing_hqr.c
 *
 * Testing file for all combinations of hierarchical QR trees.
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
        int all_llvl[] = { 0, 1, 2, 3, 4 };
        int all_hlvl[] = { 0, 1, 2, 3, 4 };
        int all_qr_p[] = { 3, 5, 7, 8 };
        int all_qr_a[] = { 1, 2, 4, 7 };
        int all_mt[]   = { 1, 3, 4, 10, 17, 25, 128 };
        int all_nt[]   = { 1, 2, 5, 13, 26, 58 };

        int nb_llvl = sizeof( all_llvl ) / sizeof( int );
        int nb_hlvl = sizeof( all_hlvl ) / sizeof( int );
        int nb_qr_p = sizeof( all_qr_p ) / sizeof( int );
        int nb_qr_a = sizeof( all_qr_a ) / sizeof( int );
        int nb_mt   = sizeof( all_mt   ) / sizeof( int );
        int nb_nt   = sizeof( all_nt   ) / sizeof( int );

        int l, h, p, a, m, n;
        int done, todo;
        todo = 0;
        done = 0;

        /* HQR */
        todo += nb_llvl * nb_mt * nb_nt * (2 * nb_qr_a - 1) * (1 +  (2 * nb_hlvl - 1) * nb_qr_p);

        for( l=0; l<nb_llvl; l++) {
            llvl = all_llvl[l];

            /*
             * Tests w/o high level trees
             */
            hlvl   = -1;
            domino = 0;
            qr_p   = -1;
            matrix.nodes = 1;
            matrix.p = 1;

            /* qr_a */
            for( a=0; a<nb_qr_a; a++) {
                qr_a = all_qr_a[a];

                /* MT */
                for( m=0; m<nb_mt; m++) {
                    MT = all_mt[m];
                    matrix.mt = MT;

                    /* NT */
                    for( n=0; n<nb_nt; n++) {
                        NT = all_nt[n];
                        matrix.nt = NT;

                        /* TS Round-robin */
                        for( tsrr=0; tsrr<2; tsrr++) {

                            /* Useless case */
                            if ((tsrr == 1) &&
                                (a    == 1))
                            {
                                continue;
                            }

                            rc = libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix,
                                                  llvl, hlvl, qr_a, qr_p, domino, tsrr );
                            if ( rc == 0 ) {
                                rc = libhqr_check( &qrtree );
                                libhqr_finalize( &qrtree );
                            }

                            if (rc != 0) {
                                fprintf(stderr, "%s -M %d -N %d --treel=%d --qr_a=%d %s      FAILED(%d)\n",
                                        argv[0], MT, NT, llvl, qr_a, (tsrr ? "-r" : "  "), rc);
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

            /*
             * Tests w/ high level trees
             */
            for( domino=0; domino<2; domino++) { /* Domino */
                if ((domino == 1) && (llvl == LIBHQR_GREEDY1P_TREE)) {
                    continue;
                }

                /* High level tree */
                for( h=0; h<nb_hlvl; h++) {
                    hlvl = all_hlvl[h];

                    /* qr_p */
                    for( p=0; p<nb_qr_p; p++) {
                        qr_p = all_qr_p[p];

                        P = qr_p;
                        matrix.nodes = P;
                        matrix.p     = P;

                        /* qr_a */
                        for( a=0; a<nb_qr_a; a++) {
                            qr_a = all_qr_a[a];

                            /* MT */
                            for( m=0; m<nb_mt; m++) {
                                MT = all_mt[m];
                                matrix.mt = MT;

                                /* NT */
                                for( n=0; n<nb_nt; n++) {
                                    NT = all_nt[n];
                                    matrix.nt = NT;

                                    /* TS Round-robin */
                                    for( tsrr=0; tsrr<2; tsrr++) {

                                        /* Useless case */
                                        if ((tsrr == 1) &&
                                            (a    == 1))
                                        {
                                            continue;
                                        }

                                        rc = libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix,
                                                         llvl, hlvl, qr_a, qr_p, domino, tsrr );

                                        if ( rc == 0 ) {
                                            rc = libhqr_check( &qrtree );
                                            libhqr_finalize( &qrtree );
                                        }

                                        if (rc != 0) {
                                            fprintf(stderr, "%s -M %d -N %d --treel=%d --treeh=%d --qr_a=%d --qr_p=%d -d=%d %s      FAILED(%d)\n",
                                                    argv[0], MT, NT, llvl, hlvl, qr_a, qr_p, domino, (tsrr ? "-r" : "  "), rc);
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
                }
            }
        }
    }
    else {
        matrix.nodes = nodes;
        matrix.p     = P;
        matrix.mt    = MT;
        matrix.nt    = NT;

        rc = libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix,
                              llvl, hlvl, qr_a, qr_p, domino, tsrr );

        if ( rc == 0 ) {
            rc = libhqr_check( &qrtree );
            libhqr_finalize( &qrtree );
        }

        if (rc != 0) {
            fprintf(stderr, "%s -M %d -N %d --treel=%d --treeh=%d --qr_a=%d --qr_p=%d -d=%d %s      FAILED(%d)\n",
                    argv[0], MT, NT, llvl, hlvl, qr_a, qr_p, domino, (tsrr ? "-r" : "  "), rc);
            ret++;
        }
        else {
            fprintf(stderr, "%s -M %d -N %d --treel=%d --treeh=%d --qr_a=%d --qr_p=%d -d=%d %s      SUCCESS\n",
                    argv[0], MT, NT, llvl, hlvl, qr_a, qr_p, domino, (tsrr ? "-r" : "  "));
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
