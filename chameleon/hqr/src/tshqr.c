/**
 *
 * @file tshqr.c
 *
 * @copyright 2018-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2018-09-30
 *
 */
#include "libhqr_internal.h"
#include "libhqr_list.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

static int
tshqr_getnbgeqrf( const libhqr_tree_t *qrtree, int k )
{
    (void)qrtree;
    (void)k;
    return 1;
}

static int
tshqr_getm( const libhqr_tree_t *qrtree, int k, int i )
{
    libhqr_tile_args_t *args = (libhqr_tile_args_t *)( qrtree->args );
    (void)i;
    return args->nbgeqrf == NULL ? 0 : k;
}

int
libhqr_initmtx_tshqr( libhqr_tree_t         *qrtree,
                      int                    roundrobin,
                      libhqr_facto_e         trans,
                      int                    kt,
                      int                    lt,
                      const libhqr_matrix_t *A )
{
    libhqr_tile_args_t *args;
    libhqr_tile_info_t *tileinfo;
    int                *pivots;
    int                 j, k, l, killer;

    if ( qrtree == NULL ) {
        fprintf( stderr, "libhqr_hqr_init, illegal value of qrtree" );
        return -1;
    }
    if ( ( trans != LIBHQR_TSQR ) && ( trans != LIBHQR_TSLQ ) ) {
        fprintf( stderr, "libhqr_hqr_init, illegal value of trans" );
        return -2;
    }
    if ( A == NULL ) {
        fprintf( stderr, "libhqr_hqr_init, illegal value of A" );
        return -3;
    }

    memset( qrtree, 0, sizeof( libhqr_tree_t ) );
    libhqr_rdmtx_initfct( qrtree );
    qrtree->init       = LIBHQR_QRTREE_MTX;
    qrtree->facto      = trans;
    qrtree->getnbgeqrf = tshqr_getnbgeqrf;
    qrtree->getm       = tshqr_getm;
    qrtree->mt         = ( trans == LIBHQR_TSQR ) ? A->mt : A->nt;
    qrtree->nt         = kt;

    qrtree->p = 1;
    qrtree->a = 1;

    qrtree->args   = malloc( sizeof( libhqr_tile_args_t ) );
    args           = (libhqr_tile_args_t *)qrtree->args;
    args->tileinfo = malloc( qrtree->mt * qrtree->nt * sizeof( libhqr_tile_info_t ) );
    args->killers  = NULL;
    args->pivots   = malloc( qrtree->nt * sizeof( int ) );
    args->nbgeqrf  = roundrobin ? malloc( sizeof(int) ) : NULL;
    tileinfo       = args->tileinfo;
    pivots         = args->pivots;

    /* Initialize the matrix */
    for ( k = 0; k < qrtree->nt; k++, tileinfo += qrtree->mt ) {
        killer = roundrobin ? k : 0;

        tileinfo[killer].type          = LIBHQR_KILLED_BY_LOCALTREE;
        tileinfo[killer].index         = killer;
        tileinfo[killer].currpiv       = -1;
        tileinfo[killer].nextpiv       = qrtree->mt;
        tileinfo[killer].prevpiv       = qrtree->mt;
        tileinfo[killer].first_nextpiv = qrtree->mt;
        tileinfo[killer].first_prevpiv = qrtree->mt;

        l = killer;
        for ( j = 1; j < qrtree->mt - lt; j++ ) {
            l = (l+1) % (qrtree->mt - lt);

            tileinfo[l].type          = LIBHQR_KILLED_BY_TS;
            tileinfo[l].index         = -2;
            tileinfo[l].currpiv       = killer;
            tileinfo[l].nextpiv       = qrtree->mt;
            tileinfo[l].prevpiv       = qrtree->mt;
            tileinfo[l].first_nextpiv = qrtree->mt;
            tileinfo[l].first_prevpiv = qrtree->mt;

            /* If first_nextpiv is undefined, it is the current tile */
            if ( tileinfo[killer].first_nextpiv == qrtree->mt ) {
                tileinfo[killer].first_nextpiv = l;
            }

            /* The previous tile is the former first_prevpiv, and the current one become the
             * first_prevpiv */
            if ( tileinfo[killer].first_prevpiv != qrtree->mt ) {
                tileinfo[tileinfo[killer].first_prevpiv].nextpiv = l;
                tileinfo[l].prevpiv                              = tileinfo[killer].first_prevpiv;
            }
            tileinfo[killer].first_prevpiv = l;
        }

        /* Add the new element when performing a TT kernel */
        if ( lt > 0 ) {
            int last = qrtree->mt - lt;

            for ( j = last; j < qrtree->mt; j++ ) {
                tileinfo[j].type          = -1;
                tileinfo[j].index         = -1;
                tileinfo[j].currpiv       = -1;
                tileinfo[j].nextpiv       = -1;
                tileinfo[j].prevpiv       = -1;
                tileinfo[j].first_nextpiv = -1;
                tileinfo[j].first_prevpiv = -1;
            }
            lt--;
        }

        /* Let's register elt1 as the first pivot */
        pivots[k] = killer;
    }

    return 0;
}

int
libhqr_init_tshqr( libhqr_tree_t         *qrtree,
                   int                    roundrobin,
                   libhqr_facto_e         trans,
                   int                    kt,
                   int                    lt,
                   const libhqr_matrix_t *A )
{
    return libhqr_initmtx_tshqr( qrtree, roundrobin, trans, kt, lt, A );
}
