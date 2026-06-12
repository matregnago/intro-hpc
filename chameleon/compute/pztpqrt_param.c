/**
 *
 * @file pztpqrt_param.c
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
#include <stdio.h>
#include "libhqr.h"

#define ATop(m,n) ATop, (m), (n)
#define A(m,n) A, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)

/**
 *  Parallel tile QR matrix reduction - Equivalent to tpqrt kernel for matrices.
 *
 * @param[in] genD
 *         Indicate if copies of the geqrt tiles must be done to speedup
 *         computations in updates. genD is considered only if D is not NULL.
 *
 * @param[in] uplo
 *         - ChamUpper: This corresponds to the former TTQRT kernel. Only the upper
 *           trapezoidal part of A is factorized.
 *         - ChamLower, or ChamUpperLower: This corresponds to the former TSQRT
 *           kernel. The full A is factorized.
 */
void chameleon_pztpqrt_param( int genD, cham_uplo_t uplo, int K,
                              const libhqr_tree_t *qrtree, CHAM_desc_t *ATop, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, n, p;
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
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        p = chameleon_pzgeqrf_param_step( genD, uplo, k, ib, qrtree, tiles,
                                          A, TS, TT, D, &options, sequence );

        /* Combine with ATop and A by merging last pivot with A(k,k) */
        {
            CHAM_desc_t *T;
            int temppm = ATop->get_blkdim( ATop, p, DIM_m, ATop->m );
            int tempkn = ATop->get_blkdim( ATop, k, DIM_n, ATop->n );
            int L, node, tempnn;

            T = TT;
            L = temppm;

            node = A->get_rankof( A, p, k );
            RUNTIME_data_migrate( sequence, ATop(k, k), node );
            RUNTIME_data_migrate( sequence, A(   p, k), node );

            INSERT_TASK_ztpqrt(
                &options,
                temppm, tempkn, chameleon_min(L, tempkn), ib, T->nb,
                ATop(k, k),
                A(p, k),
                T(p, k));

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                node = A->get_rankof( A, p, n );
                RUNTIME_data_migrate( sequence, ATop(k, n), node );
                RUNTIME_data_migrate( sequence, A(   p, n), node );

                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamConjTrans,
                    temppm, tempnn, A->nb, L, ib, T->nb,
                    A(p, k),
                    T(p, k),
                    ATop(k, n),
                    A(p, n));
            }

            chameleon_data_flush( sequence, A(p, k), request->flush );
            chameleon_data_flush( sequence, T(p, k), request->flush );
        }

        /* Restore the original location of the tiles */
        for (n = k; n < ATop->nt; n++) {
            RUNTIME_data_migrate( sequence, ATop(k, n),
                                  ATop->get_rankof( ATop, k, n ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
