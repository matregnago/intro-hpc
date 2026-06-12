/**
 *
 * @file pztpgqrt_param.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define QTop(m,n) QTop, m, n

#define A(m,n) A, m, n
#define Q(m,n) Q, m, n
#define T(m,n) T, m, n
#define D(m,n) D, m, n

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if the copies of the A tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 * @param[in] uplo
 *         - ChamUpper: This corresponds to the TTQRT factorization kernel. Only the upper
 *           trapezoidal part of Q is referenced.
 *         - ChamLower, or ChamUpperLower: This corresponds to the TSQRT factorization
 *           kernel. The full Q is referenced.
 */
void chameleon_pztpgqrt_param( int genD, cham_uplo_t uplo, int K,
                               const libhqr_tree_t *qrtree,
                               CHAM_desc_t *A, CHAM_desc_t *QTop, CHAM_desc_t *Q,
                               CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, n, p;
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
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /*
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    for (k = K-1; k >=0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        /* Setting the order of tiles */
        nbtiles = libhqr_walk_stepk( qrtree, k, tiles );
        p = tiles[nbtiles];

        /* Combine Bottom and Top matrices by merging last pivot with ATop(k,*) */
        {
            CHAM_desc_t *T = TT;
            int temppm = Q->get_blkdim( Q, p, DIM_m, Q->m );
            int tempkn = A->get_blkdim( A, k, DIM_n, A->n );
            int tempnn;
            int L = temppm;

            for (n = k; n < Q->nt; n++) {
                tempnn = Q->get_blkdim( Q, n, DIM_n, Q->n );

                node = Q->get_rankof( Q, p, n );
                RUNTIME_data_migrate( sequence, QTop(k, n), node );
                RUNTIME_data_migrate( sequence, Q(   p, n), node );

                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamNoTrans,
                    temppm, tempnn, tempkn, L, ib, T->nb,
                    A(   p, k),
                    T(   p, k),
                    QTop(k, n),
                    Q(   p, n));
            }

            chameleon_data_flush( sequence, A(p, k), request->flush );
            chameleon_data_flush( sequence, T(p, k), request->flush );
        }

        chameleon_pzungqr_param_step( genD, uplo, k, ib,
                                      qrtree, nbtiles, tiles,
                                      A, Q, TS, TT, D,
                                      &options, sequence );

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
