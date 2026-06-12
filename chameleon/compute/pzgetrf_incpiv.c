/**
 *
 * @file pzgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv parallel algorithm
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
 * @author Alycia Lisito
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
//ALLOC_WS :  ib*L->nb
//WS_ADD :  ib*L->nb
#include "control/common.h"
#include <stdlib.h>

#define A(_m_,_n_) A, _m_, _n_
#if defined(CHAMELEON_COPY_DIAG)
#define D(k)   D, k, 0
#else
#define D(k)   A, k, k
#endif
#define L(_m_,_n_) L,  _m_,  _n_
#define IPIV(_m_,_n_) &(IPIV[(int64_t)A->mb*((int64_t)(_m_)+(int64_t)A->mt*(int64_t)(_n_))])

/**
 *  Parallel tile LU factorization - dynamic scheduling
 */
void chameleon_pzgetrf_incpiv( CHAM_desc_t *A, CHAM_desc_t *L, CHAM_desc_t *D, int *IPIV,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int minMNT = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( D == NULL ) {
        D = A;
    }

    /*
     * zgetrf_incpiv = 0
     * zgessm        = 0
     * ztstrf        = A->mb * ib
     * zssssm        = 0
     */
    ws_worker = A->mb * ib;

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < minMNT; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        INSERT_TASK_zgetrf_incpiv(
            &options,
            tempkm, tempkn, ib, L->nb,
            A(k, k),
            L(k, k),
            IPIV(k, k),
            k == A->mt-1, A->nb*k);

        if ( k < (minMNT-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            INSERT_TASK_zlacpy(
                &options,
                ChamUpperLower, tempkm, tempkn,
                A(k, k),
                D(k));
#endif
        }

        for (n = k+1; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            INSERT_TASK_zgessm(
                &options,
                tempkm, tempnn, tempkm, ib, L->nb,
                IPIV(k, k),
                L(k, k),
                D(k),
                A(k, n));
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_ztstrf(
                &options,
                tempmm, tempkn, ib, L->nb,
                A(k, k),
                A(m, k),
                L(m, k),
                IPIV(m, k),
                m == A->mt-1, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_zssssm(
                    &options,
                    A->nb, tempnn, tempmm, tempnn, A->nb, ib, L->nb,
                    A(k, n),
                    A(m, n),
                    L(m, k),
                    A(m, k),
                    IPIV(m, k));
            }
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
