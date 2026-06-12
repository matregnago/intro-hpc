/**
 *
 * @file pzgeqrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the geqrt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 */
int chameleon_pzgeqrf_step( int genD, int k, int ib,
                            CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                            RUNTIME_option_t *options, RUNTIME_sequence_t *sequence )
{
    RUNTIME_request_t *request = options->request;
    int m, n;
    int tempkm, tempkn, tempnn, tempmm;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    INSERT_TASK_zgeqrt(
        options,
        tempkm, tempkn, ib, T->nb,
        A(k, k),
        T(k, k));

    if ( genD ) {
        int tempDkm = D->get_blkdim( D, k, DIM_m, D->m );
        int tempDkn = D->get_blkdim( D, k, DIM_n, D->n );
        INSERT_TASK_zlacpy(
            options,
            ChamLower, tempDkm, tempDkn,
            A(k, k),
            D(k) );
#if defined(CHAMELEON_USE_CUDA)
        INSERT_TASK_zlaset(
            options,
            ChamUpper, tempDkm, tempDkn,
            0., 1.,
            D(k) );
#endif
    }
    for (n = k+1; n < A->nt; n++) {
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        INSERT_TASK_zunmqr(
            options,
            ChamLeft, ChamConjTrans,
            tempkm, tempnn, tempkm, ib, T->nb,
            D(k),
            T(k, k),
            A(k, n));
    }
    chameleon_data_flush( sequence, D(k)   , request->flush );
    chameleon_data_flush( sequence, T(k, k), request->flush );

    for (m = k+1; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        RUNTIME_data_migrate( sequence, A(k, k),
                              A->get_rankof( A, m, k ) );

        /* TS kernel */
        INSERT_TASK_ztpqrt(
            options,
            tempmm, tempkn, 0, ib, T->nb,
            A(k, k),
            A(m, k),
            T(m, k));

        for (n = k+1; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, m, n ) );

            /* TS kernel */
            INSERT_TASK_ztpmqrt(
                options,
                ChamLeft, ChamConjTrans,
                tempmm, tempnn, A->nb, 0, ib, T->nb,
                A(m, k),
                T(m, k),
                A(k, n),
                A(m, n));
        }
        chameleon_data_flush( sequence, A(m, k), request->flush );
        chameleon_data_flush( sequence, T(m, k), request->flush );
    }

    return 1;
}

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the geqrt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 */
void chameleon_pzgeqrf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, n, ib, minMNT;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        minMNT = A->nt;
    } else {
        minMNT = A->mt;
    }

    if ( (genD == 0) || (D == NULL) ) {
        D    = A;
        genD = 0;
    }

    /*
     * zgeqrt  = A->nb * (ib+1)
     * zunmqr  = A->nb * ib
     * ztpqrt  = A->nb * (ib+1)
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < minMNT; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        chameleon_pzgeqrf_step( genD, k, ib,
                                A, T, D, &options, sequence );

        /* Restore the original location of the tiles */
        for (n = k; n < A->nt; n++) {
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
