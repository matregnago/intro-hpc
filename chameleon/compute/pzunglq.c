/**
 *
 * @file pzunglq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define Q(m,n) Q,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void chameleon_pzunglq( int genD, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempnn, tempmm, tempkmin, tempkn;
    int tempAkm, tempAkn;
    int ib, minMT;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        minMT = A->nt;
    } else {
        minMT = A->mt;
    }

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmlq  = A->nb * ib
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = minMT-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempAkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempAkn = A->get_blkdim( A, k, DIM_n, A->n );
        tempkmin = chameleon_min( tempAkn, tempAkm );
        tempkn = Q->get_blkdim( Q, k, DIM_n, Q->n );

        for (n = Q->nt-1; n > k; n--) {
            tempnn = Q->get_blkdim( Q, n, DIM_n, Q->n );
            for (m = k; m < Q->mt; m++) {
                tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

                RUNTIME_data_migrate( sequence, Q(m, k),
                                      Q->get_rankof( Q, m, n ) );

                /* TS kernel */
                INSERT_TASK_ztpmlqt(
                    &options,
                    ChamRight, ChamNoTrans,
                    tempmm, tempnn, tempAkm, 0, ib, T->nb,
                    A(k, n),
                    T(k, n),
                    Q(m, k),
                    Q(m, n));
            }
            chameleon_data_flush( sequence, A(k, n), request->flush );
            chameleon_data_flush( sequence, T(k, n), request->flush );
        }

        if ( genD ) {
            int tempDkn = D->get_blkdim( D, k, DIM_n, D->n );
            INSERT_TASK_zlacpy(
                &options,
                ChamUpper, tempkmin, tempDkn,
                A(k, k),
                D(k) );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                &options,
                ChamLower, tempkmin, tempDkn,
                0., 1.,
                D(k) );
#endif
        }
        for (m = k; m < Q->mt; m++) {
            tempmm = Q->get_blkdim( Q, m, DIM_m, Q->m );

            /* Restore the original location of the tiles */
            RUNTIME_data_migrate( sequence, Q(m, k),
                                  Q->get_rankof( Q, m, k ) );

            INSERT_TASK_zunmlq(
                &options,
                ChamRight, ChamNoTrans,
                tempmm, tempkn, tempkmin, ib, T->nb,
                D(k),
                T(k, k),
                Q(m, k));
        }
        chameleon_data_flush( sequence, D(k)   , request->flush );
        chameleon_data_flush( sequence, T(k, k), request->flush );

        RUNTIME_iteration_pop(chamctxt);
    }
    /* Mark written data for synchronization */
    Q->sync = 1;

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
