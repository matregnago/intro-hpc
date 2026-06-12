/**
 *
 * @file pzgelqf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf_param parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
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
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the gelqt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 * @param[in] uplo
 *         - ChamUpper: Classic LQ factorization of the matrix A.
 *         - ChamLower: LQ factorization of the TTLQT kernel.
 *         - ChamUpperLower: LQ factorization of the TSLQT kernel.
 */
int chameleon_pzgelqf_param_step( int genD, cham_uplo_t uplo, int k, int ib,
                                  const libhqr_tree_t *qrtree, int *tiles,
                                  CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                                  RUNTIME_option_t *options, RUNTIME_sequence_t *sequence )
{
    RUNTIME_request_t *request = options->request;
    CHAM_desc_t *T;
    int m, n, i, p;
    int L, nbgelqt;
    int tempkmin, tempkm, tempnn, tempmm, temppn;
    int node, nbtiles;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );

    /* The number of geqrt to apply */
    nbgelqt = qrtree->getnbgeqrf( qrtree, k );

    T = TS;
    for (i = 0; i < nbgelqt; i++) {
        p = qrtree->getm( qrtree, k, i );

        /* We skip the LQ factorization if this is the last diagonal tile */
        if ( (uplo == ChamLower) && (p == k) ) {
            continue;
        }

        temppn = A->get_blkdim( A, p, DIM_n, A->n );
        tempkmin = chameleon_min(tempkm, temppn);

        INSERT_TASK_zgelqt(
            options,
            tempkm, temppn, ib, T->nb,
            A(k, p), T(k, p));

        if ( genD ) {
            int tempDkm = D->get_blkdim( D, k, DIM_m, D->m );
            int tempDpn = D->get_blkdim( D, p, DIM_n, D->n );

            INSERT_TASK_zlacpy(
                options,
                ChamUpper, tempDkm, tempDpn,
                A(k, p), D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                options,
                ChamLower, tempDkm, tempDpn,
                0., 1.,
                D(k, p) );
#endif
        }

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zunmlq(
                options,
                ChamRight, ChamConjTrans,
                tempmm, temppn, tempkmin, ib, T->nb,
                D(k, p),
                T(k, p),
                A(m, p));
        }

        if ( genD || ((k+1) < A->mt)) {
            chameleon_data_flush( sequence, D(k, p), request->flush );
        }
        chameleon_data_flush( sequence, T(k, p), request->flush );
    }

    /* Setting the order of the tiles */
    nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

    for (i = 0; i < nbtiles; i++) {
        n = tiles[i];
        p = qrtree->currpiv( qrtree, k, n );

        tempnn = A->get_blkdim( A, n, DIM_n, A->n );

        if ( qrtree->gettype( qrtree, k, n ) == LIBHQR_KILLED_BY_TS ) {
            /* TS kernel */
            T = TS;
            L = 0;

            /* Force TT kernel if this is the last diagonal tile */
            if ( (uplo == ChamLower) && (n == k) ) {
                L = tempnn;
            }
        }
        else {
            /* TT kernel */
            T = TT;
            L = tempnn;
        }

        node = A->get_rankof( A, k, n );
        RUNTIME_data_migrate( sequence, A(k, p), node );
        RUNTIME_data_migrate( sequence, A(k, n), node );

        INSERT_TASK_ztplqt(
            options,
            tempkm, tempnn, chameleon_min(L, tempkm), ib, T->nb,
            A(k, p),
            A(k, n),
            T(k, n));

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            node = A->get_rankof( A, m, n );
            RUNTIME_data_migrate( sequence, A(m, p), node );
            RUNTIME_data_migrate( sequence, A(m, n), node );

            INSERT_TASK_ztpmlqt(
                options,
                ChamRight, ChamConjTrans,
                tempmm, tempnn, tempkm, L, ib, T->nb,
                A(k, n),
                T(k, n),
                A(m, p),
                A(m, n));
        }
        chameleon_data_flush( sequence, A(k, n), request->flush );
        chameleon_data_flush( sequence, T(k, n), request->flush );
    }

    return tiles[nbtiles];
}

/**
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if copies of the gelqt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 */
void chameleon_pzgelqf_param( int genD, int K,
                              const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m;
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
     * zgelqt  = A->nb * (ib+1)
     * zunmlq  = A->nb * ib
     * ztplqt  = A->nb * (ib+1)
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /*
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        chameleon_pzgelqf_param_step( genD, ChamUpper, k, ib, qrtree, tiles,
                                      A, TS, TT, D, &options, sequence );

        /* Restore the original location of the tiles */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
