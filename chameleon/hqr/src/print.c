/**
 *
 * @file print.c
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
#include "libhqr_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/**
 * @brief Print the matrix of tiles' type in the reduction trees
 * @param[in] qrtree
 *         The reduction tree to print
 */
void
libhqr_print_type( const libhqr_tree_t *qrtree )
{
    int minMN = qrtree->nt;
    int m, k;
    int lmg = 0;
    int rank = 0;

    printf("\n------------ Localization = Type of pivot --------------\n");
    for(m=0; m<qrtree->mt; m++) {
        printf("%3d | ", m);
        for (k=0; k<libhqr_imin(minMN, m+1); k++) {
            printf( "%3d ", qrtree->gettype( qrtree, k, m ) );
        }
        for (k=libhqr_imin(minMN, m+1); k<minMN; k++) {
            printf( "    " );
        }

        printf("    ");
        printf("%2d,%3d | ", rank, lmg);
        for (k=0; k<libhqr_imin(minMN, lmg+1); k++) {
            printf( "%3d ", qrtree->gettype( qrtree, k, lmg) );
        }
        for (k=libhqr_imin(minMN, lmg+1); k<minMN; k++) {
            printf( "    " );
        }
        lmg += qrtree->p;
        if ( lmg >= qrtree->mt ) {
            rank++;
            lmg = rank;
        }
        printf("\n");
    }
}

/**
 * @brief Print the pivot of each tile
 * @param[in] qrtree
 *         The reduction tree to print
 */
void
libhqr_print_pivot( const libhqr_tree_t *qrtree )
{
    int minMN = qrtree->nt;
    int m, k;
    int lmg = 0;
    int rank = 0;
    printf("\n------------ Current Pivot--------------\n");
    for(m=0; m<qrtree->mt; m++) {
        printf("%3d | ", m);
        for (k=0; k<libhqr_imin(minMN, m+1); k++) {
            printf( "%3d ", qrtree->currpiv(qrtree, k, m) );
        }
        for (k=libhqr_imin(minMN, m+1); k<minMN; k++) {
            printf( "    " );
        }

        printf("    ");
        printf("%2d,%3d | ", rank, lmg);
        for (k=0; k<libhqr_imin(minMN, lmg+1); k++) {
            printf( "%3d ", qrtree->currpiv(qrtree, k, lmg) );
        }
        for (k=libhqr_imin(minMN, lmg+1); k<minMN; k++) {
            printf( "    " );
        }
        lmg += qrtree->p;
        if ( lmg >= qrtree->mt ) {
            rank++;
            lmg = rank;
        }
        printf("\n");
    }
}

/**
 * @brief Print the next tile killed by a given pivot at step k
 * @param[in] qrtree
 *         The reduction tree to print
 */
void
libhqr_print_next_k( const libhqr_tree_t *qrtree, int k )
{
    int m, s;
    printf("\n------------ Next (k = %d)--------------\n", k);

    printf( "      " );
    for(s=qrtree->mt; s>0; s--) {
        printf( "%3d ", s );
    }
    printf( "\n" );

    for(m=0; m<qrtree->mt; m++) {
        printf("%3d | ", m);
        for(s=qrtree->mt; s>0; s--) {
            printf( "%3d ", qrtree->nextpiv(qrtree, k, m, s) );
        }
        printf("\n");
    }
}

/**
 * @brief Print the previous tile killed by a given pivot at step k
 * @param[in] qrtree
 *         The reduction tree to print
 * @parma[in] k
 *         The step k to print
 */
void
libhqr_print_prev_k( const libhqr_tree_t *qrtree, int k )
{
    int m, s;
    printf("\n------------ Prev (k = %d)--------------\n", k);

    printf( "      " );
    for(s=qrtree->mt; s>-1; s--) {
        printf( "%3d ", s );
    }
    printf( "\n" );

    for(m=0; m<qrtree->mt; m++) {
        printf("%3d | ", m);
        for(s=qrtree->mt; s>-1; s--) {
            printf( "%3d ", qrtree->prevpiv(qrtree, k, m, s) );
        }
        printf("\n");
    }
}

/**
 * @brief Print the number of geqrt per column
 * @param[in] qrtree
 *         The reduction tree to print
 * @parma[in] k
 *         The step k to print
 */
void
libhqr_print_nbgeqrt( const libhqr_tree_t *qrtree )
{
    int k, minMN = qrtree->nt;

    printf("\n------------ Nb GEQRT per k --------------\n");
    printf(" k      : ");
    for (k=0; k<minMN; k++) {
        printf( "%3d ", k );
    }
    printf( "\n" );
    printf(" Compute: ");
    for (k=0; k<minMN; k++) {
        int m, nb = 0;
        for (m=k; m < qrtree->mt; m++) {
            if ( qrtree->gettype(qrtree, k, m) > 0 ) {
                nb++;
            }
        }
        printf( "%3d ", nb );
    }
    printf( "\n" );
    printf(" Formula: ");
    for (k=0; k<minMN; k++) {
        printf( "%3d ", qrtree->getnbgeqrf( qrtree, k ) );
    }
    printf( "\n" );
}

/**
 * @brief Print the list of geqrt indices at step k
 * @param[in] qrtree
 *         The reduction tree to print
 * @parma[in] k
 *         The step k to print
 */
void
libhqr_print_geqrt_k( const libhqr_tree_t *qrtree, int k )
{
    int i, nb;

    printf("\n------------ Liste of geqrt for k = %d --------------\n", k);

    printf( "  m:");
    nb = qrtree->getnbgeqrf( qrtree, k );
    for (i=0; i < nb; i++) {
        int m = qrtree->getm( qrtree, k, i );
        if ( i == qrtree->geti( qrtree, k, m) ) {
            printf( "%3d ", m );
        }
        else {
            printf( "x%2d ", qrtree->geti( qrtree, k, m) );
        }
    }
    printf( "\n" );
}
