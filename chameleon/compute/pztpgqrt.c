/**
 *
 * @file pztpgqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define V2(m,n) V2,  m,  n
#define T2(m,n) T2,  m,  n
#define Q1(m,n) Q1,  m,  n
#define Q2(m,n) Q2,  m,  n
#define D(k)    D,   k,  k

/**
 *  Parallel tile QR factorization - dynamic scheduling
 */
void chameleon_pztpgqrt( int KT, int L,
                         CHAM_desc_t *V2, CHAM_desc_t *T2,
                         CHAM_desc_t *Q1, CHAM_desc_t *Q2,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempkn, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = chameleon_max( Q2->m - L, 1 );
    int maxmt = (maxm % Q2->mb == 0) ? (maxm / Q2->mb) : (maxm / Q2->mb + 1);
    int maxmtk;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    /*
     * ztpmqrt = Q1->nb * ib
     */
    ws_worker = Q1->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * ztpmqrt = 3 * Q1->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * Q1->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = KT-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkn = Q1->get_blkdim( Q1, k, DIM_n, Q1->n );

        /* Equivalent to the tsmqr step on Q1,Q2 */
        maxmtk = chameleon_min( Q2->mt, maxmt+k ) - 1;
        for (m = maxmtk; m > -1; m--) {
            tempmm = Q2->get_blkdim( Q2, m, DIM_m, Q2->m );
            templm = ((L > 0) && (m == maxmtk)) ? tempmm : 0;

            for (n = k; n < Q2->nt; n++) {
                tempnn = Q2->get_blkdim( Q2, n, DIM_n, Q2->n );
                /* TT kernel */
                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamNoTrans,
                    tempmm, tempnn, tempkn, templm, ib, T2->nb,
                    V2(m, k),
                    T2(m, k),
                    Q1(k, n),
                    Q2(m, n) );
            }
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
