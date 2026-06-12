/**
 *
 * @file pzgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv parallel algorithm
 *
 * @version 1.4.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Terry Cojean
 * @author Matthieu Kuhn
 * @author Pierre Esterie
 * @author Alycia Lisito
 * @author Xavier Lacoste
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n)  A,  m, n
#define WD(m)    WL, m, m
#define WL(m, n) WL, m, n
#define WU(m, n) WU, m, n

/**
 * @brief Generic tile algorithm of the LU factorization without pivoting
 *
 * This is the version to use by default.
 */
void chameleon_pzgetrf_nopiv_generic( CHAM_desc_t        *A,
                                      RUNTIME_sequence_t *sequence,
                                      RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n, ib;
    int tempkm, tempkn, tempmm, tempnn;
    int min_mnt = chameleon_min( A->mt, A->nt );
    int kmin, kmax;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

#if defined(CHAMELEON_USE_MPI)
    /*
     * Estimate the number of tasks per step on each node to automatically limit
     * the submission window and prevent the memory overflow issue dur to
     * pre-allocation of the reception buffers.
     */
    if ( chamctxt->autominmax_enabled && (chamctxt->scheduler == RUNTIME_SCHED_STARPU) ) {
        int P                = chameleon_desc_datadist_get_iparam(A, 0);
        int Q                = chameleon_desc_datadist_get_iparam(A, 1);
        int lookahead        = chamctxt->lookahead;
        int nbtasks_per_step = (A->mt * A->nt) / (P * Q);
        int mintasks         = nbtasks_per_step *  lookahead;
        int maxtasks         = nbtasks_per_step * (lookahead+1);

        if ( CHAMELEON_Comm_rank() == 0 ) {
            chameleon_warning( "chameleon_pzgetrf_nopiv",
                               "Setting limit for the number of submitted tasks\n" );
        }
        RUNTIME_set_minmax_submitted_tasks( mintasks, maxtasks );
    }
#endif

    kmin = chameleon_max( 0,       chamctxt->first_step );
    kmax = chameleon_min( min_mnt, chamctxt->last_step  );
    for ( k = kmin; k < kmax; k++ ) {

        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );

        options.priority = 2*A->nt - 2*k;
        INSERT_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), A->mb*k);

        for ( m = k+1; m < A->mt; m++ ) {
            options.priority = 2*A->nt - 2*k - m;
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_ztrsm(
                &options,
                ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
                tempmm, tempkn, A->mb,
                zone, A(k, k),
                      A(m, k));
        }
        for ( n = k+1; n < A->nt; n++ ) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            options.priority = 2*A->nt - 2*k - n;
            INSERT_TASK_ztrsm(
                &options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, A(k, k),
                      A(k, n));

            for ( m = k+1; m < A->mt; m++ ) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                options.priority = 2*A->nt - 2*k  - n - m;
                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, A(m, k),
                           A(k, n),
                    zone,  A(m, n));
            }

            chameleon_data_flush( sequence, A(k, n), request->flush );
        }

        for (m = k; m < A->mt; m++) {
            chameleon_data_flush( sequence, A(m, k), request->flush );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_finalize(&options, chamctxt);
}

/**
 * @brief Tile algorithm of the LU factorization without pivoting using
 * workspace to optimize the communications
 *
 * This version should be used only the workspaces have been initialized. It
 * used a ring of communication to propagate the column and row panel at each
 * iteration to regulate the flow of tasks.
 * By doing so, the row and column panel are communicated along a ring with a
 * given lookahead. Thus the number of step of the algorithm ongoing at given
 * instant `t` is limited to lookahead steps.
 */
void chameleon_pzgetrf_nopiv_ws( CHAM_desc_t        *A,
                                 CHAM_desc_t        *WL,
                                 CHAM_desc_t        *WU,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int k, m, n, ib, lp, lq;
    int tempkm, tempkn, tempmm, tempnn;
    int lookahead, myp, myq, P, Q;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;
    lookahead = chamctxt->lookahead;
    P   = chameleon_desc_datadist_get_iparam(A, 0);
    Q   = chameleon_desc_datadist_get_iparam(A, 1);
    myp = A->myrank / Q;
    myq = A->myrank % Q;

    for ( k = 0; k < chameleon_min( A->mt, A->nt ); k++ ) {
        RUNTIME_iteration_push(chamctxt, k);
        lp = (k % lookahead) * P;
        lq = (k % lookahead) * Q;

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );

        options.priority = 2*A->nt - 2*k;
        INSERT_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), A->mb*k);

        /**
         * Broadcast of A(k,k) along rings in both directions
         */
        chameleon_pzbcast_tile( ChamRowwise, ChamBcastRing,
                                A, k, k, WL, k, lq, &options );
        chameleon_pzbcast_tile( ChamColumnwise, ChamBcastRing,
                                A, k, k, WU, lp, k, &options );

        chameleon_data_flush( sequence, A( k, k ), request->flush );

        for ( m = k+1; m < A->mt; m++ ) {

            /* Skip the row if you are not involved with */
            if ( ( m % P ) != myp ) {
                continue;
            }

            options.priority = 2*A->nt - 2*k - m;
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            assert( A->get_rankof( A, m, k ) == WU->get_rankof( WU, myp + lp, k) );
            INSERT_TASK_ztrsm(
                &options,
                ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
                tempmm, tempkn, A->mb,
                zone, WU( myp + lp, k ),
                      A(  m,        k ) );

            /* Broadcast A(m,k) into temp buffers through a ring */
            chameleon_pzbcast_tile( ChamRowwise, ChamBcastRing,
                                    A, m, k, WL, m, lq, &options );

            chameleon_data_flush( sequence, A( m, k ), request->flush );
        }

        for ( n = k+1; n < A->nt; n++ ) {

            /* Skip the column if you are not involved with */
            if ( ( n % Q ) != myq ) {
                continue;
            }

            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            options.priority = 2*A->nt - 2*k - n;

            assert( A->get_rankof( A, k, n ) == WL->get_rankof( WL, k, myq+lq) );
            INSERT_TASK_ztrsm(
                &options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, WL( k, myq + lq ),
                      A( k, n ));

            /* Broadcast A(k,n) into temp buffers through a ring */
            chameleon_pzbcast_tile( ChamColumnwise, ChamBcastRing,
                                    A, k, n, WU, lp, n, &options );

            chameleon_data_flush( sequence, A( k, n ), request->flush );

            for ( m = k+1; m < A->mt; m++ ) {

                /* Skip the row if you are not involved with */
                if ( ( m % P ) != myp ) {
                    continue;
                }

                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                options.priority = 2*A->nt - 2*k  - n - m;

                assert( A->get_rankof( A, m, n ) == WL->get_rankof( WL, m, myq + lq) );
                assert( A->get_rankof( A, m, n ) == WU->get_rankof( WU, myp + lp, n) );

                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, WL( m, myq + lq ),
                           WU( myp + lp, n ),
                    zone,  A( m, n ));
            }
        }
        RUNTIME_iteration_pop( chamctxt );
    }

    CHAMELEON_Desc_Flush( WL, sequence );
    CHAMELEON_Desc_Flush( WU, sequence );

    RUNTIME_options_finalize( &options, chamctxt );
}

void chameleon_pzgetrf_nopiv( struct chameleon_pzgetrf_nopiv_s *ws,
                              CHAM_desc_t                      *A,
                              RUNTIME_sequence_t               *sequence,
                              RUNTIME_request_t                *request )
{
    if ( ws && ws->use_workspace ) {
        chameleon_pzgetrf_nopiv_ws( A, &(ws->WL), &(ws->WU), sequence, request );
    }
    else {
        chameleon_pzgetrf_nopiv_generic( A, sequence, request );
    }
    /* Mark written data for synchronization */
    A->sync = 1;
}
