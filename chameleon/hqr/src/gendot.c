/**
 *
 * @file gendot.c
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

char *color[] = {
    "red", "blue", "green", "orange", "cyan", "purple", "yellow",
};
#define DAG_NBCOLORS 7

void
libhqr_print_dot( const libhqr_tree_t *qrtree, const char *filename )
{
    int  *pos, *next, *done;
    int   k, m, n, lpos, prev, length;
    int   minMN = libhqr_imin( qrtree->mt, qrtree->nt );
    FILE *f     = fopen( filename, "w" );

    if ( f == NULL ) {
        perror("fopen");
        return;
    }

    done = (int *)malloc( qrtree->mt * sizeof( int ) );
    pos  = (int *)malloc( qrtree->mt * sizeof( int ) );
    next = (int *)malloc( qrtree->mt * sizeof( int ) );
    memset( pos, 0, qrtree->mt * sizeof( int ) );
    memset( next, 0, qrtree->mt * sizeof( int ) );

    /* Print header */
    fprintf( f, "digraph G { orientation=portrait; \n" );

    /* Print labels */
    for ( m = 0; m < qrtree->mt; m++ ) {
        fprintf( f, "%d [label=\"%d\",color=white,pos=\"-1.,-%d.!\"]\n", m, m, m );
    }

    for ( k = 0; k < minMN; k++ ) {
        int nb2reduce = qrtree->mt - k - 1;

        for ( m = k; m < qrtree->mt; m++ ) {
            fprintf( f,
                     "p%d_m%d_k%d [shape=point,width=0.1, pos=\"%d.,-%d.!\",color=\"%s\"];\n",
                     m, qrtree->mt, k, pos[m], m,
                     color[( m % qrtree->p ) % DAG_NBCOLORS] );
            next[m] = qrtree->nextpiv( qrtree, k, m, qrtree->mt );
        }

        while ( nb2reduce > 0 ) {
            memset( done, 0, qrtree->mt * sizeof( int ) );
            for ( m = qrtree->mt - 1; m > ( k - 1 ); m-- ) {
                n = next[m];
                if ( next[n] != qrtree->mt ) {
                    continue;
                }
                if ( n != qrtree->mt ) {
                    lpos = libhqr_imax( pos[m], pos[n] );
                    lpos++;
                    pos[m] = lpos;
                    pos[n] = lpos;

                    fprintf( f, "p%d_m%d_k%d [shape=point,width=0.1, pos=\"%d.,-%d.!\",color=\"%s\"];\n",
                             m, n, k, pos[m], m,
                             color[( m % qrtree->p ) % DAG_NBCOLORS] );

                    prev = qrtree->prevpiv( qrtree, k, m, n );
                    fprintf( f, "p%d_m%d_k%d->p%d_m%d_k%d [width=0.1,color=\"%s\"]\n",
                             m, prev, k, m, n, k,
                             color[( m % qrtree->p ) % DAG_NBCOLORS] );

                    prev = qrtree->prevpiv( qrtree, k, n, n );
                    if ( qrtree->gettype( qrtree, k, n ) == 0 ) {
                        fprintf( f, "p%d_m%d_k%d->p%d_m%d_k%d [style=dotted, width=0.1,color=\"%s\"]\n",
                                 n, prev, k, m, n, k,
                                 color[( m % qrtree->p ) % DAG_NBCOLORS] );
                    }
                    else {
                        fprintf( f, "p%d_m%d_k%d->p%d_m%d_k%d [style=dashed, width=0.1,color=\"%s\"]\n",
                                 n, prev, k, m, n, k,
                                 color[( m % qrtree->p ) % DAG_NBCOLORS] );
                    }

                    next[m] = qrtree->nextpiv( qrtree, k, m, n );
                    done[m] = done[n] = 1;
                    nb2reduce--;
                }
            }
        }
    }

    length = 0;
    for ( m = 0; m < qrtree->mt; m++ ) {
        length = libhqr_imax( length, pos[m] );
    }
    length++;
    for ( k = 0; k < length; k++ ) {
        fprintf( f, "l%d [label=\"%d\",color=white,pos=\"%d.,0.5!\"]\n", k, k, k );
    }
    fprintf( f, "} // close graph\n" );

    printf( "Tic Max = %d\n", length - 1 );

    fclose( f );
    free( done );
    free( pos );
    free( next );
}
