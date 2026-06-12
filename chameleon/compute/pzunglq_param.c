/**
 *
 * @file pzunglq_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq_param parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m,n) A, (m), (n)
#define Q(m,n) Q, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)

/**
 *  Parallel construction of Q using tile V - dynamic scheduling
 */
void chameleon_pzunglq_param( int genD, const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int K, L;
    int tempkm, tempkmin, temppn, tempnn, tempmm;
    int ib, node, nbtiles, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmlq  = A->nb * ib
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int));

    K = chameleon_min(A->mt, A->nt);

    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );

        /* Setting the order of the tiles*/
        nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

        for (i = nbtiles-1; i >= 0; i--) {
            n = tiles[i];
            p = qrtree->currpiv(qrtree, k, n);

            tempnn = Q->get_blkdim( Q, n, DIM_n, Q->n );

            if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                /* TS kernel */
                L = 0;
                T = TS;
            }
            else {
                /* TT kernel */
                L = tempnn;
                T = TT;
            }
            for (m = k; m < Q->mt; m++) {
                tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                node = Q->get_rankof( Q, m, n );
                RUNTIME_data_migrate( sequence, Q(m, p), node );
                RUNTIME_data_migrate( sequence, Q(m, n), node );

                INSERT_TASK_ztpmlqt(
                    &options,
                    ChamRight, ChamNoTrans,
                    tempmm, tempnn, tempkm, L, ib, T->nb,
                    A(k, n),
                    T(k, n),
                    Q(m, p),
                    Q(m, n));
            }
            chameleon_data_flush( sequence, A(k, n), request->flush );
            chameleon_data_flush( sequence, T(k, n), request->flush );
        }

        T = TS;
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            p = qrtree->getm(qrtree, k, i);

            temppn = A->get_blkdim( A, p, DIM_n, A->n );
            tempkmin = chameleon_min(tempkm, temppn);

            if ( genD ) {
                int tempDpn = D->get_blkdim( D, p, DIM_n, D->n );
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpper, tempkmin, tempDpn,
                    A(k, p),
                    D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamLower, tempkmin, tempDpn,
                    0., 1.,
                    D(k, p) );
#endif
            }
            for (m = k; m < Q->mt; m++) {
                tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                RUNTIME_data_migrate( sequence, Q(m, p),
                                      Q->get_rankof( Q, m, p ) );

                INSERT_TASK_zunmlq(
                    &options,
                    ChamRight, ChamNoTrans,
                    tempmm, temppn, tempkmin, ib, T->nb,
                    D(k, p),
                    T(k, p),
                    Q(m, p));
            }
            chameleon_data_flush( sequence, D(k, p), request->flush );
            chameleon_data_flush( sequence, T(k, p), request->flush );
        }

        RUNTIME_iteration_pop(chamctxt);
    }
    /* Mark written data for synchronization */
    Q->sync = 1;

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
