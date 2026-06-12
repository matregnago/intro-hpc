/**
 *
 * @file tphqr.c
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

#if defined(LIBHQR_DEBUG)
static inline void
print_onelist( libhqr_list_t *list )
{
    libhqr_list_elt_t *elt = *list;
    while( elt != NULL ) {
        printf("{ %-3d, %-3d } ", elt->id, elt->date );
        elt = elt->next;
    }
    printf("\n");
}

static inline void
print_alllist( int nblist, libhqr_list_t *lists )
{
    int i;

    for(i=0; i<nblist; i++) {
        print_onelist( lists + i );
    }
}
#endif

static int
tphqr_getnbgeqrf( const libhqr_tree_t *qrtree, int k )
{
    libhqr_tile_args_t *args    = (libhqr_tile_args_t*)(qrtree->args);
    int                *nbgeqrf = args->nbgeqrf;
    return nbgeqrf[k];
}

static int
tphqr_getm( const libhqr_tree_t *qrtree, int k, int i )
{
    libhqr_tile_args_t *args    = (libhqr_tile_args_t*)(qrtree->args);
    int                *killers = args->killers;

    killers += k * qrtree->mt;

    return killers[i];
}

int
libhqr_initmtx_tphqr( libhqr_tree_t         *qrtree,
                      libhqr_facto_e         trans,
                      int                    kt,
                      int                    lt,
                      const libhqr_matrix_t *A,
                      libhqr_tree_e          hlvl,
                      int                    a,
                      int                    p )
{
    libhqr_tile_args_t *args;
    libhqr_tile_info_t *tileinfo;
    int                *killers, *pivots;
    int                 low_mt, i, j, k;
    int                 lp, la;

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

    /* Compute parameters */
    if ( a == -1 ) {
        a = 4;
    }
    else {
        a = libhqr_imax( a, 1 );
    }

    if ( p == -1 ) {
        if ( trans == LIBHQR_TSQR ) {
            p = A->p;
        }
        else {
            p = A->nodes / A->p;
        }
    }
    else {
        p = libhqr_imax( p, 1 );
    }

    memset( qrtree, 0, sizeof( libhqr_tree_t ) );
    libhqr_rdmtx_initfct( qrtree );
    qrtree->init       = LIBHQR_QRTREE_MTX;
    qrtree->facto      = trans;
    qrtree->getnbgeqrf = tphqr_getnbgeqrf;
    qrtree->getm       = tphqr_getm;
    qrtree->mt         = ( trans == LIBHQR_TSQR ) ? A->mt : A->nt;
    qrtree->nt         = kt;

    a = libhqr_imin( a, qrtree->mt );

    low_mt = libhqr_iceil( qrtree->mt, p );
    low_mt = libhqr_iceil( low_mt, a );

    qrtree->p = p;
    qrtree->a = a;

    qrtree->args   = malloc( sizeof( libhqr_tile_args_t ) );
    args           = (libhqr_tile_args_t *)qrtree->args;
    args->tileinfo = malloc( qrtree->mt * qrtree->nt * sizeof( libhqr_tile_info_t ) );
    args->killers  = malloc( qrtree->mt * qrtree->nt * sizeof( int ) );
    args->pivots   = malloc( 2 * qrtree->nt * sizeof( int ) );
    args->nbgeqrf  = args->pivots + qrtree->nt;
    tileinfo       = args->tileinfo;
    killers        = args->killers;
    pivots         = args->pivots;

    /**
     * Compute the number of sets:
     *   -1) No list, Rows have been eliminated
     *    0) 1 for the top level (Reduction between nodes)
     *    1) p for the middle level (reduction within a node of the masters of the TS domains)
     *    2) low_mt * p for the low level TS trees
     *    3) low_mt * p to temporary save the node that will be reduced at the next step
     */
    libhqr_list_t     *lists = calloc( ( 2 * low_mt + 1 ) * p + 1, sizeof( libhqr_list_t ) );
    libhqr_list_t     *lists0, *lists1, *lists2, *lists3;
    libhqr_list_elt_t *elts = malloc( qrtree->mt * sizeof( libhqr_list_elt_t ) );
    libhqr_list_elt_t *elt1, *elt2;

    lists3 = lists;
    lists2 = lists3 + low_mt * p;
    lists1 = lists2 + low_mt * p;
    lists0 = lists1 + p;

    for ( k = 0; k < qrtree->mt - lt; k++ ) {
        elts[k].next = NULL;
        elts[k].id   = k;
        elts[k].date = 0;

        lp = k % p;
        la = ( k / p ) / a;
        assert( lp * low_mt + la < low_mt * p );
        libhqr_list_push( lists3 + lp * low_mt + la, elts + k );
    }

    /* Initialize the matrix */
    for ( k = 0; k < qrtree->nt; k++, killers += qrtree->mt, tileinfo += qrtree->mt ) {
        int l, killer_index = 0;

        /* Move remaining rows from level 3 to level 2 */
        for ( i = 0; i < ( low_mt * p ); i++ ) {
            assert( lists2[i] == NULL );
            lists2[i] = lists3[i];
            lists3[i] = NULL;
        }

        /* Handle TS tiles */
        l = 0;
        for ( i = 0; i < p; i++ ) {
            for ( j = 0; j < low_mt; j++, l++ ) {
                elt1 = libhqr_list_pop( lists2 + l );
                elt2 = libhqr_list_pop( lists2 + l );

                if ( elt1 == NULL ) {
                    continue;
                }

                assert( ( elt1 != NULL ) && ( elt1->next == NULL ) );

                tileinfo[elt1->id].type          = LIBHQR_KILLED_BY_LOCALTREE;
                tileinfo[elt1->id].index         = killer_index;
                tileinfo[elt1->id].currpiv       = -1;
                tileinfo[elt1->id].nextpiv       = qrtree->mt;
                tileinfo[elt1->id].prevpiv       = qrtree->mt;
                tileinfo[elt1->id].first_nextpiv = qrtree->mt;
                tileinfo[elt1->id].first_prevpiv = qrtree->mt;
                killers[killer_index]            = elt1->id;
                killer_index++;

                /* Let's kill the tiles with the main one */
                while ( elt2 != NULL ) {
                    assert( ( elt2 != NULL ) && ( elt2->next == NULL ) );

                    int date   = libhqr_imax( elt1->date, elt2->date ) + 1;
                    elt1->date = date;
                    elt2->date = date;

                    tileinfo[elt2->id].type          = LIBHQR_KILLED_BY_TS;
                    tileinfo[elt2->id].index         = -2;
                    tileinfo[elt2->id].currpiv       = elt1->id;
                    tileinfo[elt2->id].nextpiv       = qrtree->mt;
                    tileinfo[elt2->id].prevpiv       = qrtree->mt;
                    tileinfo[elt2->id].first_nextpiv = qrtree->mt;
                    tileinfo[elt2->id].first_prevpiv = qrtree->mt;

                    /* If first_nextpiv is undefined, it is the current tile */
                    if ( tileinfo[elt1->id].first_nextpiv == qrtree->mt ) {
                        tileinfo[elt1->id].first_nextpiv = elt2->id;
                    }

                    /* The previous tile is the former first_prevpiv, and the current one become the
                     * first_prevpiv */
                    if ( tileinfo[elt1->id].first_prevpiv != qrtree->mt ) {
                        tileinfo[tileinfo[elt1->id].first_prevpiv].nextpiv = elt2->id;
                        tileinfo[elt2->id].prevpiv = tileinfo[elt1->id].first_prevpiv;
                    }
                    tileinfo[elt1->id].first_prevpiv = elt2->id;

                    /* Push back elt2 to level3 for next round */
                    libhqr_list_push( lists3 + l, elt2 );
                    elt2 = libhqr_list_pop( lists2 + l );
                }

                libhqr_list_push( lists1 + i, elt1 );
            }
        }

        /* Handle local tree */
        for ( i = 0; i < p; i++ ) {
            elt1 = libhqr_list_pop( lists1 + i );
            elt2 = libhqr_list_pop( lists1 + i );

            if ( elt1 == NULL ) {
                continue;
            }

            while ( elt2 != NULL ) {
                int date   = libhqr_imax( elt1->date, elt2->date ) + 1;
                elt1->date = date;
                elt2->date = date;

                assert( tileinfo[elt2->id].type == LIBHQR_KILLED_BY_LOCALTREE );
                tileinfo[elt2->id].currpiv = elt1->id;

                /* If first_nextpiv is undefined, it is the current tile */
                if ( tileinfo[elt1->id].first_nextpiv == qrtree->mt ) {
                    tileinfo[elt1->id].first_nextpiv = elt2->id;
                }

                /* The previous tile is the former first_prevpiv, and the current one become the
                 * first_prevpiv */
                if ( tileinfo[elt1->id].first_prevpiv != qrtree->mt ) {
                    tileinfo[tileinfo[elt1->id].first_prevpiv].nextpiv = elt2->id;
                    tileinfo[elt2->id].prevpiv = tileinfo[elt1->id].first_prevpiv;
                }
                tileinfo[elt1->id].first_prevpiv = elt2->id;

                /* Push back elt2 to level3 for next round */
                l = ( elt2->id / p ) / a;
                libhqr_list_push( lists3 + i * low_mt + l, elt2 );

                /* Push back elt1 */
                libhqr_list_push( lists1 + i, elt1 );
                elt1 = libhqr_list_pop( lists1 + i );
                elt2 = libhqr_list_pop( lists1 + i );
            }

            /* Upgrade last element */
            tileinfo[elt1->id].type = LIBHQR_KILLED_BY_DISTTREE;
            libhqr_list_push( lists0, elt1 );
        }

        /* Handle distributed tree */
        {
            elt1 = libhqr_list_pop( lists0 );
            elt2 = libhqr_list_pop( lists0 );

            if ( elt1 == NULL ) {
                continue;
            }

            while ( elt2 != NULL ) {
                int date   = libhqr_imax( elt1->date, elt2->date ) + 1;
                elt1->date = date;
                elt2->date = date;

                assert( tileinfo[elt2->id].type == LIBHQR_KILLED_BY_DISTTREE );
                tileinfo[elt2->id].currpiv = elt1->id;

                /* If first_nextpiv is undefined, it is the current tile */
                if ( tileinfo[elt1->id].first_nextpiv == qrtree->mt ) {
                    tileinfo[elt1->id].first_nextpiv = elt2->id;
                }

                /* The previous tile is the former first_prevpiv, and the current one become the
                 * first_prevpiv */
                if ( tileinfo[elt1->id].first_prevpiv != qrtree->mt ) {
                    tileinfo[tileinfo[elt1->id].first_prevpiv].nextpiv = elt2->id;
                    tileinfo[elt2->id].prevpiv = tileinfo[elt1->id].first_prevpiv;
                }
                tileinfo[elt1->id].first_prevpiv = elt2->id;

                /* Push back elt2 to level3 for next round */
                i = elt2->id % p;
                l = ( elt2->id / p ) / a;
                libhqr_list_push( lists3 + i * low_mt + l, elt2 );

                if ( hlvl != LIBHQR_FLAT_TREE ) {
                    libhqr_list_push( lists0, elt1 );
                    elt1 = libhqr_list_pop( lists0 );
                }
                elt2 = libhqr_list_pop( lists0 );
            }

            i = elt1->id % p;
            l = ( elt1->id / p ) / a;
            libhqr_list_push( lists3 + i * low_mt + l, elt1 );

            /* Let's register elt1 as the first pivot */
            pivots[k] = elt1->id;
        }

        assert( killer_index <= qrtree->mt );
        args->nbgeqrf[k] = killer_index;

        /* Add the new element when performing a TT kernel */
        if ( lt > 0 ) {
            int last        = qrtree->mt - lt;
            elts[last].next = NULL;
            elts[last].id   = last;
            elts[last].date = 0;

            lp = last % p;
            la = ( last / p ) / a;
            assert( lp * low_mt + la < low_mt * p );
            libhqr_list_push( lists3 + lp * low_mt + la, elts + last );

            for ( i = last; i < qrtree->mt; i++ ) {
                tileinfo[i].type          = -1;
                tileinfo[i].index         = -1;
                tileinfo[i].currpiv       = -1;
                tileinfo[i].nextpiv       = -1;
                tileinfo[i].prevpiv       = -1;
                tileinfo[i].first_nextpiv = -1;
                tileinfo[i].first_prevpiv = -1;
            }
            lt--;
        }
    }

#if defined(LIBHQR_DEBUG)
    if ( 0 ) {
        tileinfo = args->tileinfo;
        for ( j = 0; j < qrtree->nt; j++, tileinfo += qrtree->mt ) {
            printf( "---- Step %d ----\n", j );
            for ( i = 0; i < qrtree->mt; i++ ) {
                printf( "{type = %-1d, index = %-2d, currpiv = %-3d, nextpiv = %-3d, prevpiv = "
                        "%-3d, first_nextpiv = %-3d, first_prevpiv = %-3d}\n",
                        tileinfo[i].type,
                        tileinfo[i].index,
                        tileinfo[i].currpiv,
                        tileinfo[i].nextpiv,
                        tileinfo[i].prevpiv,
                        tileinfo[i].first_nextpiv,
                        tileinfo[i].first_prevpiv );
            }
        }
    }

    if ( 0 ) {
        print_alllist( low_mt * p, lists );
    }
#endif

    free( elts );
    free( lists );
    return 0;
}

int
libhqr_init_tphqr( libhqr_tree_t         *qrtree,
                   libhqr_facto_e         trans,
                   int                    kt,
                   int                    lt,
                   const libhqr_matrix_t *A,
                   libhqr_tree_e          hlvl,
                   int                    a,
                   int                    p )
{
    return libhqr_initmtx_tphqr( qrtree, trans, kt, lt, A, hlvl, a, p );
}
