/**
 *
 * @file pzunglqrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglqrh parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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

#define A(m,n)  A,  (m),  (n)
#define Q(m,n)  Q,  (m),  (n)
#define T(m,n)  T,  (m),  (n)
#define T2(m,n) T,  (m),  (n)+(A->nt)
#define D(m,n)  D,  (m),  (n)

/**
 *  Parallel construction of Q using tile V (application to identity;
 *  reduction Householder) - dynamic scheduling
 */
void chameleon_pzunglqrh( int genD, int BS,
                          CHAM_desc_t *A, CHAM_desc_t *Q,
                          CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, N, RD, lastRD;
    int tempkm, tempkmin, tempNn, tempnn, tempmm, tempNRDn;
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

    K = chameleon_min(A->mt, A->nt);
    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        lastRD = 0;
        for (RD = BS; RD < A->nt-k; RD *= 2)
            lastRD = RD;
        for (RD = lastRD; RD >= BS; RD /= 2) {
            for (N = k; N+RD < A->nt; N += 2*RD) {
                tempNRDn = A->get_blkdim( A, N+RD, DIM_n, A->n );
                for (m = k; m < Q->mt; m++) {
                    tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                    node = Q->get_rankof( Q, m, N+RD );
                    RUNTIME_data_migrate( sequence, Q(m, N),    node );
                    RUNTIME_data_migrate( sequence, Q(m, N+RD), node );

                    /* TT kernel */
                    INSERT_TASK_ztpmlqt(
                        &options,
                        ChamRight, ChamNoTrans,
                        tempmm, tempNRDn, tempkm, tempNRDn, ib, T->nb,
                        A (k, N+RD),
                        T2(k, N+RD),
                        Q (m, N   ),
                        Q (m, N+RD));
                }

                chameleon_data_flush( sequence, A (k, N+RD), request->flush );
                chameleon_data_flush( sequence, T2(k, N+RD), request->flush );
            }
        }
        for (N = k; N < A->nt; N += BS) {
            tempNn = A->get_blkdim( A, N, DIM_n, A->n );
            tempkmin = chameleon_min(tempkm, tempNn);
            for (n = chameleon_min(N+BS, A->nt)-1; n > N; n--) {
                tempnn = Q->get_blkdim( Q, n, DIM_n, Q->n );

                for (m = k; m < Q->mt; m++) {
                    tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                    node = Q->get_rankof( Q, m, n );
                    RUNTIME_data_migrate( sequence, Q(m, N), node );
                    RUNTIME_data_migrate( sequence, Q(m, n), node );

                    /* TS kernel */
                    INSERT_TASK_ztpmlqt(
                        &options,
                        ChamRight, ChamNoTrans,
                        tempmm, tempnn, tempkm, 0, ib, T->nb,
                        A(k, n),
                        T(k, n),
                        Q(m, N),
                        Q(m, n));
                }

                chameleon_data_flush( sequence, A(k, n), request->flush );
                chameleon_data_flush( sequence, T(k, n), request->flush );
            }

            if ( genD ) {
                int tempDNn = D->get_blkdim( D, N, DIM_n, D->n );

                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpper, tempkmin, tempDNn,
                    A(k, N),
                    D(k, N) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamLower, tempkmin, tempDNn,
                    0., 1.,
                    D(k, N) );
#endif
            }
            for (m = k; m < Q->mt; m++) {
                tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                RUNTIME_data_migrate( sequence, Q(m, N),
                                      Q->get_rankof( Q, m, N ) );

                INSERT_TASK_zunmlq(
                    &options,
                    ChamRight, ChamNoTrans,
                    tempmm, tempNn,
                    tempkmin, ib, T->nb,
                    D(k, N),
                    T(k, N),
                    Q(m, N));
            }
            chameleon_data_flush( sequence, D(k, N), request->flush );
            chameleon_data_flush( sequence, T(k, N), request->flush );
        }
        RUNTIME_iteration_pop(chamctxt);
    }
    /* Mark written data for synchronization */
    Q->sync = 1;

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
