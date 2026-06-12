/**
 *
 * @file pzgelqfrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqfrh parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
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
#define T(m,n)  T,  (m),  (n)
#define T2(m,n) T,  (m),  (n)+A->nt
#define D(m,n)  D,  (m),  (n)

/*
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 */
void chameleon_pzgelqfrh( int genD, int BS, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, N, RD;
    int tempkmin, tempkm, tempNn, tempnn, tempmm, tempNRDn;
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
     * zgelqt  = A->nb * (ib+1)
     * zunmlq  = A->nb * ib
     * ztplqt  = A->nb * (ib+1)
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
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

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );

        for (N = k; N < A->nt; N += BS) {
            tempNn = A->get_blkdim( A, N, DIM_n, A->n );
            tempkmin = chameleon_min(tempkm, tempNn);
            INSERT_TASK_zgelqt(
                &options,
                tempkm, tempNn, ib, T->nb,
                A(k, N),
                T(k, N));
            if ( genD ) {
                int tempDkm = D->get_blkdim( D, k, DIM_m, D->m );
                int tempDNn = D->get_blkdim( D, N, DIM_n, D->n );

                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpper, tempDkm, tempDNn,
                    A(k, N),
                    D(k, N) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamLower, tempDkm, tempDNn,
                    0., 1.,
                    D(k, N) );
#endif
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_zunmlq(
                    &options,
                    ChamRight, ChamConjTrans,
                    tempmm, tempNn, tempkmin, ib, T->nb,
                    D(k, N),
                    T(k, N),
                    A(m, N));
            }
            chameleon_data_flush( sequence, D(k, N), request->flush );
            chameleon_data_flush( sequence, T(k, N), request->flush );

            for (n = N+1; n < chameleon_min(N+BS, A->nt); n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                RUNTIME_data_migrate( sequence, A(k, N),
                                      A->get_rankof( A, k, n ) );

                /* TS kernel */
                INSERT_TASK_ztplqt(
                    &options,
                    tempkm, tempnn, 0, ib, T->nb,
                    A(k, N),
                    A(k, n),
                    T(k, n));

                for (m = k+1; m < A->mt; m++) {
                    tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                    RUNTIME_data_migrate( sequence, A(m, N),
                                          A->get_rankof( A, m, n ) );

                    INSERT_TASK_ztpmlqt(
                        &options,
                        ChamRight, ChamConjTrans,
                        tempmm, tempnn, tempkm, 0, ib, T->nb,
                        A(k, n),
                        T(k, n),
                        A(m, N),
                        A(m, n));
                }
                chameleon_data_flush( sequence, A(k, n), request->flush );
                chameleon_data_flush( sequence, T(k, n), request->flush );
            }
        }
        for (RD = BS; RD < A->nt-k; RD *= 2) {
            for (N = k; N+RD < A->nt; N += 2*RD) {
                tempNRDn = A->get_blkdim( A, N+RD, DIM_n, A->n );

                node = A->get_rankof( A, k, N+RD );
                RUNTIME_data_migrate( sequence, A(k, N),    node );
                RUNTIME_data_migrate( sequence, A(k, N+RD), node );

                /* TT kernel */
                INSERT_TASK_ztplqt(
                    &options,
                    tempkm, tempNRDn, chameleon_min(tempNRDn, tempkm), ib, T->nb,
                    A (k, N   ),
                    A (k, N+RD),
                    T2(k, N+RD));

                for (m = k+1; m < A->mt; m++) {
                    tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                    node = A->get_rankof( A, m, N+RD );
                    RUNTIME_data_migrate( sequence, A(m, N),    node );
                    RUNTIME_data_migrate( sequence, A(m, N+RD), node );

                    INSERT_TASK_ztpmlqt(
                        &options,
                        ChamRight, ChamConjTrans,
                        tempmm, tempNRDn, tempkm, tempNRDn, ib, T->nb,
                        A (k, N+RD),
                        T2(k, N+RD),
                        A (m, N   ),
                        A (m, N+RD));
                }
                chameleon_data_flush( sequence, A (k, N+RD), request->flush );
                chameleon_data_flush( sequence, T2(k, N+RD), request->flush );
            }
        }

        /* Restore the original location of the tiles */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
