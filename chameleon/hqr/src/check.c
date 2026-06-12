/**
 *
 * @file check.c
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

#define ENDCHECK( test, ret )                   \
    do {                                        \
        if ( !test ) {                          \
            assert( ret == 0 );                 \
            return ret;                         \
        }                                       \
    } while (0)

/**
 * @brief Check Formula for NB geqrt
 * @param[in] qrtree
 *         The reduction tree on which to check the number of geqrt calls.
 * @retval True if ok, False if an error was found.
 *
 * This function checks the number of tiles of type 1, 2 or 3 per column and
 * make sure it matches the value returned by getnbgeqrf() subfunction.
 */
int
libhqr_check_nbgeqrf( const libhqr_tree_t *qrtree )
{
    int check = 1;
    int k, m;
    int checkTS = (qrtree->facto == LIBHQR_TSQR) || (qrtree->facto == LIBHQR_TSLQ);

#if defined(DEBUG_WITH_PRINT)
    libhqr_print_type( qrtree );
    libhqr_print_nbgeqrt( qrtree );
#endif

    for (k=0; k<qrtree->nt; k++) {
        int nb = 0;
        int m0 = checkTS ? 0 : k;
        for (m=m0; m < qrtree->mt; m++) {
            if ( qrtree->gettype( qrtree, k, m ) > 0 ) {
                nb++;
            }
        }

        if ( nb != qrtree->getnbgeqrf( qrtree, k ) ) {
            check = 0;
            printf(" ----------------------------------------------------\n"
                   "  - a = %d, p = %d, M = %d, N = %d\n"
                   "     Check number of geqrt:\n"
                   "       For k=%d => return %d instead of %d\n",
                   qrtree->a, qrtree->p, qrtree->mt, qrtree->nt, k, qrtree->getnbgeqrf( qrtree, k ), nb );
        }
    }

    return check;
}

/**
 * @brief Check geqrt indices
 * @param[in] qrtree
 *         The reduction tree on which to check the geqrt indices.
 * @retval True if ok, False if an error was found.
 *
 * This function checks the bijection between the getm() and geti()
 * subfunctions.
 */
int
libhqr_check_geti_getm( const libhqr_tree_t *qrtree )
{
    int check = 1;
    int k, i, m;
    int checkTS = (qrtree->facto == LIBHQR_TSQR) || (qrtree->facto == LIBHQR_TSLQ);

    for (k=0; k<qrtree->nt; k++) {
        int nb = qrtree->getnbgeqrf( qrtree, k );
        int prevm = -1;
        for (i=0; i < nb; i++) {

            m = qrtree->getm( qrtree, k, i );

            /*
             * getm has to be the inverse of geti
             */
            if ( i != qrtree->geti( qrtree, k, m) ) {
                check = 0;
                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Check indices of geqrt:\n"
                       "        getm( k=%d, i=%d ) => m = %d && geti( k=%d, m=%d ) => i = %d\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                       k, i, m, k, m, qrtree->geti( qrtree, k, m));
            }
            /*
             * Tiles before the diagonal are factorized and
             * the m is a growing list (not true with round-robin inside TS)
             */
            else if ( !checkTS && (qrtree->a == 1) && (( m < k ) || ( m < prevm )) ) {
                check = 0;
                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Check indices of geqrt:\n"
                       "        getm( k=%d, i=%d ) => m = %d\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt, k, i, m);
            }
            prevm = m;
        }
    }

    return check;
}

/**
 * @brief Check the number of exit with next subfunction
 * @param[in] qrtree
 *         The reduction tree on which to check.
 * @retval True if ok, False if an error was found.
 *
 * This function checks that one and only one exit entry exists per column in
 * the next subfunction.
 */
int
libhqr_check_next_exit( const libhqr_tree_t *qrtree )
{
    int k, m, s;
    int check = 1;
    int checkTS = (qrtree->facto == LIBHQR_TSQR) || (qrtree->facto == LIBHQR_TSLQ);

    for (k=0; k<qrtree->nt; k++) {
        int m0 = checkTS ? 0 : k;
        for(m=m0; m<qrtree->mt; m++) {
            int nb = 0;
            for(s=qrtree->mt; s>=m0; s--) {
                if ( qrtree->nextpiv(qrtree, k, m, s) == qrtree->mt ) {
                    nb++;
                }
            }
            if ( nb > 1 ) {
                libhqr_print_next_k( qrtree, k);
                libhqr_print_prev_k( qrtree, k);

                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Next of line %d for step %d contains more than one exit:\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                       m, k);
                check = 0;
            }
            else if ( nb == 0 ) {
                libhqr_print_next_k( qrtree, k);
                libhqr_print_prev_k( qrtree, k);

                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Next of line %d for step %d needs one exit:\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                       m, k);
                check = 0;
            }
        }
    }
    return check;
}

/**
 * @brief Check the number of exit with prev subfunction
 * @param[in] qrtree
 *         The reduction tree on which to check.
 * @retval True if ok, False if an error was found.
 *
 * This function checks that one and only one exit entry exists per column in
 * the prev subfunction.
 */
int
libhqr_check_prev_exit( const libhqr_tree_t *qrtree )
{
    int k, m, s;
    int check = 1;
    int checkTS = (qrtree->facto == LIBHQR_TSQR) || (qrtree->facto == LIBHQR_TSLQ);

    for (k=0; k<qrtree->nt; k++) {

        int m0 = checkTS ? 0 : k;
        for(m=m0; m<qrtree->mt; m++) {
            int nb = 0;
            for(s=m0; s<qrtree->mt; s++) {
                if ( qrtree->prevpiv(qrtree, k, m, s) == qrtree->mt ) {
                    nb++;
                }
            }
            if ( nb > 1 ) {
                libhqr_print_next_k( qrtree, k);
                libhqr_print_prev_k( qrtree, k);

                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Prev of line %d for step %d contains more than one exit:\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                       m, k);
                check = 0;
            }
            else if ( nb == 0 ) {
                libhqr_print_next_k( qrtree, k);
                libhqr_print_prev_k( qrtree, k);

                printf(" ----------------------------------------------------\n"
                       "  - a = %d, p = %d, M = %d, N = %d\n"
                       "     Prev of line %d for step %d needs one exit:\n",
                       qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                       m, k);
                check = 0;
            }
        }
    }
    return check;
}

/**
 * @brief Check the bijection next(prev()) = id()
 * @param[in] qrtree
 *         The reduction tree on which to check.
 * @retval True if ok, False if an error was found.
 *
 * This function checks that the subfunctions prev() and next() are bijectives.
 */
int
libhqr_check_prev_next( const libhqr_tree_t *qrtree )
{
    int k, m;
    int check = 1;
    int checkTS = (qrtree->facto == LIBHQR_TSQR) || (qrtree->facto == LIBHQR_TSLQ);

    for (k=0; k<qrtree->nt; k++) {
        int start = qrtree->mt;
        int m0 = checkTS ? 0 : k;
        int next, prev;

        for(m=m0; m<qrtree->mt; m++) {

            do {
                next = qrtree->nextpiv(qrtree, k, m, start);
                if ( next == qrtree->mt ) {
                    prev = qrtree->prevpiv(qrtree, k, m, m);
                }
                else {
                    prev = qrtree->prevpiv(qrtree, k, m, next);
                }

                if ( start != prev ) {
                    libhqr_print_next_k( qrtree, k);
                    libhqr_print_prev_k( qrtree, k);

                    printf(" ----------------------------------------------------\n"
                           "  - a = %d, p = %d, M = %d, N = %d\n"
                           "     Check next/prev:\n"
                           "       next( k=%d, m=%d, start=%d ) => %d && prev( k=%d, m=%d, start=%d ) => %d (instead of %d)\n",
                           qrtree->a, qrtree->p, qrtree->mt, qrtree->nt,
                           k, m, start, next, k, m, next, prev, start);
                    check = 0;
                }
                start = next;
            } while ( start != qrtree->mt );
        }
    }
    return check;
}

/**
 * @brief Check the correctness of a given reduction tree
 * @param[in] qrtree
 *         The reduction tree to check.
 * @retval 0 if ok, >0 if an error occured.
 */
int
libhqr_check( const libhqr_tree_t *qrtree )
{
    int check;

    /* Check the number of geqrt */
    check = libhqr_check_nbgeqrf( qrtree );
    ENDCHECK( check, 1 );

    /* Check indices of geqrt */
    check = libhqr_check_geti_getm( qrtree );
    ENDCHECK( check, 2 );

    /* Check number of exit in next */
    check = libhqr_check_next_exit( qrtree );
    ENDCHECK( check, 3 );

    /* Check number of exit in prev */
    check = libhqr_check_prev_exit( qrtree );
    ENDCHECK( check, 4 );

    /* Check next/prev */
    check = libhqr_check_prev_next( qrtree );
    ENDCHECK( check, 5 );

    return 0;
}
