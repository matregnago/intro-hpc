/**
 *
 * @file pzgeqrf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf_param parallel algorithm
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
#include "libhqr.h"

#define A(m,n) A, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the geqrt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 * @param[in] uplo
 *         - ChamLower: Classic QR factorization of the matrix A.
 *         - ChamUpper: QR factorization of the TTQRT kernel.
 *         - ChamUpperLower: QR factorization of the TSQRT kernel.
 */
int chameleon_pzgeqrf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                  const libhqr_tree_t *qrtree, int *tiles,
                                  CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                  RUNTIME_option_t *options, RUNTIME_sequence_t *sequence )
{
    RUNTIME_request_t *request = options->request;
    CHAM_desc_t *T;
    int m, n, i, p;
    int L, nbgeqrt;
    int tempkmin, tempkn, tempnn, tempmm, temppm;
    int node, nbtiles;

    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* The number of geqrt to apply */
    nbgeqrt = qrtree->getnbgeqrf( qrtree, k );

    T = TS;
    for (i = 0; i < nbgeqrt; i++) {
        p = qrtree->getm( qrtree, k, i );

        /* We skip the QR factorization if this is the last diagonal tile */
        if ( (uplo == ChamUpper) && (p == k) ) {
            continue;
        }

        temppm = A->get_blkdim( A, p, DIM_m, A->m );
        tempkmin = chameleon_min(temppm, tempkn);

        INSERT_TASK_zgeqrt(
            options,
            temppm, tempkn, ib, T->nb,
            A(p, k), T(p, k) );

        if ( genD ) {
            int tempDpm = D->get_blkdim( D, p, DIM_m, D->m );
            int tempDkn = D->get_blkdim( D, k, DIM_n, D->n );

            INSERT_TASK_zlacpy(
                options,
                ChamLower, tempDpm, tempDkn,
                A(p, k), D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                options,
                ChamUpper, tempDpm, tempDkn,
                0., 1.,
                D(p, k) );
#endif
        }

        for (n = k+1; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            INSERT_TASK_zunmqr(
                options,
                ChamLeft, ChamConjTrans,
                temppm, tempnn, tempkmin, ib, T->nb,
                D(p, k),
                T(p, k),
                A(p, n));
        }

        if ( genD || ((k+1) < A->nt)) {
            chameleon_data_flush( sequence, D(p, k), request->flush );
        }
        chameleon_data_flush( sequence, T(p, k), request->flush );
    }

    /* Setting the order of the tiles */
    nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

    for (i = 0; i < nbtiles; i++) {
        m = tiles[i];
        p = qrtree->currpiv( qrtree, k, m );

        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        if ( qrtree->gettype( qrtree, k, m ) == LIBHQR_KILLED_BY_TS ) {
            /* TS kernel */
            T = TS;
            L = 0;

            /* Force TT kernel if this is the last diagonal tile */
            if ( (uplo == ChamUpper) && (m == k) ) {
                L = tempmm;
            }
        }
        else {
            /* TT kernel */
            T = TT;
            L = tempmm;
        }

        node = A->get_rankof( A, m, k );
        RUNTIME_data_migrate( sequence, A(p, k), node );
        RUNTIME_data_migrate( sequence, A(m, k), node );

        INSERT_TASK_ztpqrt(
            options,
            tempmm, tempkn, chameleon_min(L, tempkn), ib, T->nb,
            A(p, k),
            A(m, k),
            T(m, k));

        for (n = k+1; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            node = A->get_rankof( A, m, n );
            RUNTIME_data_migrate( sequence, A(p, n), node );
            RUNTIME_data_migrate( sequence, A(m, n), node );

            INSERT_TASK_ztpmqrt(
                options,
                ChamLeft, ChamConjTrans,
                tempmm, tempnn, tempkn, L, ib, T->nb,
                A(m, k),
                T(m, k),
                A(p, n),
                A(m, n));
        }
        chameleon_data_flush( sequence, A(m, k), request->flush );
        chameleon_data_flush( sequence, T(m, k), request->flush );
    }

    return tiles[nbtiles];
}

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the geqrt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 */
void chameleon_pzgeqrf_param( int genD, int K,
                              const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, n;
    int ib, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

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
    /*
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push( chamctxt, k );

        chameleon_pzgeqrf_param_step( genD, ChamLower, k, ib, qrtree, tiles,
                                      A, TS, TT, D, &options, sequence );

        /* Restore the original location of the tiles */
        for (n = k; n < A->nt; n++) {
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );
        }

        RUNTIME_iteration_pop( chamctxt );
    }

    free( tiles );
    RUNTIME_options_ws_free( &options );
    RUNTIME_options_finalize( &options, chamctxt );
}
