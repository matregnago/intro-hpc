/**
 *
 * @file pztrsmpl.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsmpl parallel algorithm
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
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define L(m,n) L,  m,  n
#define IPIV(m,n) &(IPIV[(int64_t)A->nb*((int64_t)(m)+(int64_t)A->mt*(int64_t)(n))])
/**
 *  Parallel forward substitution for tile LU - dynamic scheduling
 */
void chameleon_pztrsmpl( CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib, K, DIM_k;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    if ( A->m <= A->n ) {
        K     = A->m;
        DIM_k = DIM_m;
    }
    else {
        K     = A->n;
        DIM_k = DIM_n;
    }

    ib = CHAMELEON_IB;
    for (k = 0; k < chameleon_min(A->mt, A->nt); k++) {
        tempkm   = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn   = A->get_blkdim( A, k, DIM_n, A->n );
        tempkmin = A->get_blkdim( A, k, DIM_k, K    );

        for (n = 0; n < B->nt; n++) {
            tempnn = B->get_blkdim( B, n, DIM_n, B->n );
            INSERT_TASK_zgessm(
                &options,
                tempkm, tempnn, tempkmin, ib, L->nb,
                IPIV(k, k),
                L(k, k),
                A(k, k),
                B(k, n));
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            for (n = 0; n < B->nt; n++) {
                tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                INSERT_TASK_zssssm(
                    &options,
                    A->nb, tempnn, tempmm, tempnn, tempkn, ib, L->nb,
                    B(k, n),
                    B(m, n),
                    L(m, k),
                    A(m, k),
                    IPIV(m, k));
            }
        }
    }
    /* Mark written data for synchronization */
    B->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
