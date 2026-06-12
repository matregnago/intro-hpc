/**
 *
 * @file pzunmqrrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqrrh parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A, (m), (n)
#define C(m,n)  C, (m), (n)
#define T(m,n)  T, (m), (n)
#define T2(m,n) T, (m), ((n)+A->nt)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel application of Q using tile V - QR factorization (reduction
 *  Householder) - dynamic scheduling
 */
void chameleon_pzunmqrrh( int genD, int BS, cham_side_t side, cham_trans_t trans,
                          CHAM_desc_t *A, CHAM_desc_t *C, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, p;
    int KT, RD, lastRD;
    int temppm, temppn, tempkn, tempnn, tempmm, tempkmin;
    int ib, node;

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
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    KT = chameleon_min( A->mt, A->nt );
    if (side == ChamLeft ) {
        if (trans == ChamConjTrans) {
            /*
             *  ChamLeft / ChamConjTrans
             */
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = A->get_blkdim( A, k, DIM_n, A->n );

                for (p = k; p < C->mt; p += BS) {

                    temppm = C->get_blkdim( C, p, DIM_m, C->m );
                    tempkmin = chameleon_min( temppm, tempkn );


                    if ( genD ) {
                        int tempDpm = D->get_blkdim( D, p, DIM_m, D->m );

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );
                        INSERT_TASK_zunmqr(
                            &options,
                            side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(p, n));
                    }
                    chameleon_data_flush( sequence, D(p, k), request->flush );
                    chameleon_data_flush( sequence, T(p, k), request->flush );

                    for (m = p+1; m < chameleon_min(p+BS, C->mt); m++) {
                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                        for (n = 0; n < C->nt; n++) {
                            tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(m, k),
                                T(m, k),
                                C(p, n),
                                C(m, n));
                        }
                        chameleon_data_flush( sequence, A(m, k), request->flush );
                        chameleon_data_flush( sequence, T(m, k), request->flush );
                    }
                }
                for (RD = BS; RD < C->mt-k; RD *= 2) {
                    for (p = k; p+RD < C->mt; p += 2*RD) {
                        m = p+RD;

                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                        for (n = 0; n < C->nt; n++) {
                            tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, tempmm, ib, T->nb,
                                A (m, k),
                                T2(m, k),
                                C (p, n),
                                C (m, n));
                        }
                        chameleon_data_flush( sequence, A (m, k), request->flush );
                        chameleon_data_flush( sequence, T2(m, k), request->flush );
                    }
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < C->nt; n++) {
                    RUNTIME_data_migrate( sequence, C(k, n),
                                          C->get_rankof( C, k, n ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamLeft / ChamNoTrans
         */
        else {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                lastRD = 0;
                for (RD = BS; RD < C->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (p = k; p+RD < C->mt; p += 2*RD) {
                        m = p+RD;

                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                        for (n = 0; n < C->nt; n++) {
                            tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, tempmm, ib, T->nb,
                                A (m, k),
                                T2(m, k),
                                C (p, n),
                                C (m, n));
                        }
                        chameleon_data_flush( sequence, A (m, k), request->flush );
                        chameleon_data_flush( sequence, T2(m, k), request->flush );
                    }
                }
                for (p = k; p < C->mt; p += BS) {

                    for (m = chameleon_min(p+BS, C->mt)-1; m > p; m--) {
                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                        for (n = 0; n < C->nt; n++) {
                            tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(m, k),
                                T(m, k),
                                C(p, n),
                                C(m, n));
                        }
                        chameleon_data_flush( sequence, A(m, k), request->flush );
                        chameleon_data_flush( sequence, T(m, k), request->flush );
                    }

                    temppm = C->get_blkdim( C, p, DIM_m, C->m );
                    tempkmin = chameleon_min( temppm, tempkn );

                    if ( genD ) {
                        int tempDpm = D->get_blkdim( D, p, DIM_m, D->m );

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (n = 0; n < C->nt; n++) {
                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                        RUNTIME_data_migrate( sequence, C(p, n),
                                              C->get_rankof( C, p, n ) );

                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(p, n));
                    }
                    chameleon_data_flush( sequence, D(p, k), request->flush );
                    chameleon_data_flush( sequence, T(p, k), request->flush );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /*
     *  ChamRight / ChamConjTrans
     */
    else {
        if (trans == ChamConjTrans) {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = A->get_blkdim( A, k, DIM_n, A->n );

                lastRD = 0;
                for (RD = BS; RD < C->nt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (p = k; p+RD < C->nt; p += 2*RD) {
                        n = p+RD;

                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                        for (m = 0; m < C->mt; m++) {
                            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, tempmm, ib, T->nb,
                                A (n, k),
                                T2(n, k),
                                C (m, p),
                                C (m, n));
                        }
                        chameleon_data_flush( sequence, A (n, k), request->flush );
                        chameleon_data_flush( sequence, T2(n, k), request->flush );
                    }
                }
                for (p = k; p < C->nt; p += BS) {

                    for (n = chameleon_min(p+BS, C->nt)-1; n > p; n--) {

                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                        for (m = 0; m < C->mt; m++) {
                            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(n, k),
                                T(n, k),
                                C(m, p),
                                C(m, n));
                        }
                        chameleon_data_flush( sequence, A(n, k), request->flush );
                        chameleon_data_flush( sequence, T(n, k), request->flush );
                    }

                    temppn = C->get_blkdim( C, p, DIM_n, C->n );
                    tempkmin = chameleon_min( temppn, tempkn );

                    if ( genD ) {
                        int tempDpm = D->get_blkdim( D, p, DIM_m, D->m );

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                        RUNTIME_data_migrate( sequence, C(m, p),
                                              C->get_rankof( C, m, p ) );

                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(m, p));
                    }
                    chameleon_data_flush( sequence, D(p, k), request->flush );
                    chameleon_data_flush( sequence, T(p, k), request->flush );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamNoTrans
         */
        else {
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = A->get_blkdim( A, k, DIM_n, A->n );

                for (p = k; p < C->nt; p += BS) {
                    temppn = C->get_blkdim( C, p, DIM_n, C->n );
                    tempkmin = chameleon_min( temppn, tempkn );


                    if ( genD ) {
                        int tempDpm = D->get_blkdim( D, p, DIM_m, D->m );

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(m, p));
                    }
                    chameleon_data_flush( sequence, D(p, k), request->flush );
                    chameleon_data_flush( sequence, T(p, k), request->flush );

                    for (n = p+1; n < chameleon_min(p+BS,  C->nt); n++) {
                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );
                        for (m = 0; m < C->mt; m++) {
                            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, 0, ib, T->nb,
                                A(n, k),
                                T(n, k),
                                C(m, p),
                                C(m, n));
                        }
                        chameleon_data_flush( sequence, A(n, k), request->flush );
                        chameleon_data_flush( sequence, T(n, k), request->flush );
                    }
                }
                for (RD = BS; RD < C->nt-k; RD *= 2) {
                    for (p = k; p+RD < C->nt; p += 2*RD) {
                        n = p + RD;
                        tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                        for (m = 0; m < C->mt; m++) {
                            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmqrt(
                                &options, side, trans,
                                tempmm, tempnn, tempkn, tempmm, ib, T->nb,
                                A (n, k),
                                T2(n, k),
                                C (m, p),
                                C (m, n));
                        }
                        chameleon_data_flush( sequence, A (n, k), request->flush );
                        chameleon_data_flush( sequence, T2(n, k), request->flush );
                    }
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < C->mt; m++) {
                    RUNTIME_data_migrate( sequence, C(m, k),
                                          C->get_rankof( C, m, k ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /* Mark written data for synchronization */
    C->sync = 1;

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
