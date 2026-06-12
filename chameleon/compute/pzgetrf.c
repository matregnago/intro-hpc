/**
 *
 * @file pzgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf parallel algorithm
 *
 * @version 1.4.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @author Samuel Thibault
 * @author Xavier Lacoste
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include "chameleon/flops.h"

#define A(m,n)   A,                m, n
#define Up(m)    ws->Up,           m, 0
#define Wu(m,n)  ws->laswp->Wu,    m, n
#define Wl(m,n)  ws->Wl,           m, n
#define Ws(m,n)  &(ws->laswp->ws), m, n

/*
 * Function to compute the optimal batch size when generating the graph that
 * covers all the possibilities of pivoting within the factorized panel.
 *
 * @param[in] ws
 *      The workspace data structure associated to the algorithm that holds all
 *      extra information that may be needed for LU factorization
 *
 * @param[in] mt
 *      The number of tiles to factorize within the panel
 *
 * @param[in] mb
 *      The default tile size.
 *
 * @param[in] j
 *      The index of the column to factorize
 *
 * @return The optimal batch size.
 *
 */
static inline int
chameleon_pzgetrf_batch_size( const struct chameleon_pzgetrf_s *ws,
                              int mt, int mb, int j )
{
    if ( ws->batch_adaptive == 0 ) {
        return ( ( j % ws->ib ) != 0 ) ? ws->batch_size_blas2 : ws->batch_size_blas3;
    }

    double flops     = flops_zgetrf_blocked_offdiag( mb, mb, j, ws->ib );
    int    batch_max = chameleon_min( CHAMELEON_BATCH_SIZE, mt );
    int    batch_sze = batch_max;

    /**
     * First solution. (Alycia Lisito)
     *
     * This solution aims at maximizing the batch size.
     */
    if ( j != 0 ) {
        batch_sze = chameleon_min( chameleon_max( ws->flops_min / flops, 1 ), batch_max );
    }

    if ( mt % batch_sze != 0 ) {
        batch_sze = chameleon_min( chameleon_ceil( mt, mt / batch_sze ), batch_max );
    }

    /**
     * Second solution. (Mathieu Faverge)
     *
     * This solution aims at balancing the load rather than increasing the task size.
     */
    /*
     batch_sze = chameleon_ceil( mt,
                                 chameleon_max( chameleon_ceil( mt, batch_max ),
                                                ( mt * flops ) / ws->flops_min ) );
     */
    return batch_sze;
}

/*
 * All the functions below are panel factorization variant.
 * The parameters are:
 *   @param[inout] ws
 *      The data structure associated to the algorithm that holds all extra
 *      information that may be needed for LU factorization
 *
 *   @param[inout] A
 *      The descriptor of the full matrix A (not just the panel)
 *
 *   @param[inout] ipiv
 *      The descriptor of the pivot array associated to A.
 *
 *   @param[in] k
 *      The index of the column to factorize
 *
 *   @param[inout] options
 *      The runtime options data structure to pass through all insert_task calls.
 */
static inline void
chameleon_pzgetrf_panel_facto_percol( struct chameleon_pzgetrf_s *ws,
                                      CHAM_desc_t                *A,
                                      CHAM_ipiv_t                *ipiv,
                                      CHAM_desc_pivot_t          *pivot,
                                      int                         k,
                                      RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;

    /*
     * Algorithm per column with pivoting
     */
    for (h=0; h<=minmn; h++){
        INSERT_TASK_zgetrf_percol_diag(
            options,
            tempkm, tempkn, h, k * A->mb,
            A(k, k),
            ipiv, pivot );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_percol_offdiag(
                options,
                tempmm, tempkn, h, m * A->mb,
                A(m, k),
                pivot );
        }

        /* Invalidate prevpiv because we don't need anymore */
        RUNTIME_pivot_invalidate( pivot, A->myrank, h-1 );

#if defined(CHAMELEON_USE_MPI)
        /* Reduce globally (between MPI processes) */
        INSERT_TASK_zpivot_allreduce( options, A, pivot, k, h, tempkn, ws->laswp );
#endif
    }

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, A->mt - k, ipiv, k );
    RUNTIME_pivot_flushone( options->sequence, pivot, A->myrank );

    (void)ws;
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_percol_batched( struct chameleon_pzgetrf_s *ws,
                                              CHAM_desc_t                *A,
                                              CHAM_ipiv_t                *ipiv,
                                              CHAM_desc_pivot_t          *pivot,
                                              int                         k,
                                              RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;
    void **clargs = malloc( sizeof(char *) );
    memset( clargs, 0, sizeof(char *) );

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;

    /*
     * Algorithm per column with pivoting (no recursion)
     */
    /* Iterate on current panel column */
    /* Since index h scales column h-1, we need to iterate up to minmn (included) */
    for ( h = 0; h <= minmn; h++ ) {

        INSERT_TASK_zgetrf_percol_diag( options, tempkm, tempkn, h, k * A->mb, A(k, k), ipiv, pivot );

        for ( m = k+1; m < A->mt; m++ ) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_panel_offdiag_batched( options, tempmm, tempkn, h, m * A->mb,
                                                      (void *)ws, A(m, k), clargs, pivot );
        }
        INSERT_TASK_zgetrf_panel_offdiag_batched_flush( options, A, k, clargs, pivot );

        /* Invalidate prevpiv because we don't need anymore */
        RUNTIME_pivot_invalidate( pivot, A->myrank, h-1 );

#if defined(CHAMELEON_USE_MPI)
        /* Reduce globally (between MPI processes) */
        INSERT_TASK_zpivot_allreduce( options, A, pivot, k, h, tempkn, ws->laswp );
#endif
    }

    free( clargs );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, A->mt - k, ipiv, k );
    RUNTIME_pivot_flushone( options->sequence, pivot, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto_blocked( struct chameleon_pzgetrf_s *ws,
                                       CHAM_desc_t                *A,
                                       CHAM_ipiv_t                *ipiv,
                                       CHAM_desc_pivot_t          *pivot,
                                       int                         k,
                                       RUNTIME_option_t           *options )
{
    const RUNTIME_request_t *request = options->request;
    int m, h, b, nbblock, ib, last_ib, trsmUp;
    int tempkm, tempkn, tempmm, minmn;
    int rankAkk, rankAmk, myrank;
    int hmax, j;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;
    ib       = ws->ib;
    nbblock  = chameleon_ceil( minmn, ib );
    rankAkk  = A->get_rankof( A, k, k );
    myrank   = A->myrank;

    /*
     * Algorithm per column with pivoting
     */
    for ( b = 0; b < nbblock; b++ ) {
        last_ib = ( b == ( nbblock - 1 ) );
        hmax    = ( last_ib ) ? minmn + 1 - b * ib : ib;

        for ( h = 0; h < hmax; h++ ) {
            j = h + b * ib;
            trsmUp = ( ( h == 0 ) && ( j != 0 ) && ( j != minmn ) );

            if ( trsmUp ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ib, tempkn, j, ib,
                    Up(myrank),
                    pivot );
            }

            INSERT_TASK_zgetrf_blocked_diag(
                options,
                tempkm, tempkn, j, k * A->mb, ib, trsmUp,
                A(k, k), Up(rankAkk),
                ipiv, pivot );

            for ( m = k+1; m < A->mt; m++ ) {
                tempmm  = A->get_blkdim( A, m, DIM_m, A->m );
                rankAmk = A->get_rankof( A, m, k );

                INSERT_TASK_zgetrf_blocked_offdiag(
                    options,
                    tempmm, tempkn, j, m * A->mb, ib, trsmUp,
                    A(m, k), Up(rankAmk),
                    pivot );
            }

            assert( j <= minmn );

            /* Invalidate prevpiv because we don't need anymore */
            RUNTIME_pivot_invalidate( pivot, A->myrank, j-1 );

#if defined(CHAMELEON_USE_MPI)
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zpivot_allreduce( options, A, pivot, k, j, tempkn, ws->laswp );
#endif

            /*
             * Copy the row of the pivot in Up for the replicated TRSM.
             * If b is the last inner block then we don't have a TRSM.
             * If j+1 is a multiple of ib then we do the copy of the last
             * pivot row in the TRSM task.
             */
            if ( ( !last_ib ) && ( ( ( j + 1 ) % ib ) != 0 ) ) {
                INSERT_TASK_zgetrf_cpy_pivrow_in_Up( options, ws->Up, A->myrank,
                                                     k, j, ib, tempkn, pivot );
            }
        }
    }
    chameleon_data_flush( options->sequence, Up(myrank), request->flush );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, A->mt - k, ipiv, k );
    RUNTIME_pivot_flushone( options->sequence, pivot, A->myrank );
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_blocked_batched( struct chameleon_pzgetrf_s *ws,
                                               CHAM_desc_t                *A,
                                               CHAM_ipiv_t                *ipiv,
                                               CHAM_desc_pivot_t          *pivot,
                                               int                         k,
                                               RUNTIME_option_t           *options )
{
    const RUNTIME_request_t *request = options->request;
    int m, h, b, nbblock, ib, hmax, j, last_ib, trsmUp;
    int tempkm, tempkn, tempmm, minmn;
    int rankAmk, myrank;
    void **clargs = malloc( sizeof(char *) );
    memset( clargs, 0, sizeof(char *) );

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;
    ib       = ws->ib;
    nbblock  = chameleon_ceil( minmn, ib );
    myrank   = A->myrank;

    /*
     * Algorithm per column with pivoting (no recursion)
     */
    /* Iterate on current panel column */
    /* Since index h scales column h-1, we need to iterate up to minmn (included) */
    for ( b = 0; b < nbblock; b++ ) {
        last_ib = ( b == ( nbblock - 1 ) );
        hmax    = ( last_ib ) ? minmn + 1 - b * ib : ib;

        for ( h = 0; h < hmax; h++ ) {
            j =  h + b * ib;
            trsmUp = ( ( h == 0 ) && ( j != 0 ) && ( j != minmn ) );

            if ( trsmUp ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ib, tempkn, j, ib,
                    Up(myrank),
                    pivot );
            }

            ws->batch_size = chameleon_pzgetrf_batch_size( ws, A->mt - k, A->nb, j );
            for ( m = k; m < A->mt; m++ ) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                rankAmk = A->get_rankof( A, m, k );
                INSERT_TASK_zgetrf_panel_blocked_batched( options, tempmm, tempkn, j, m * A->mb, trsmUp,
                                                          (void *)ws, A(m, k), Up(rankAmk), clargs, ipiv, pivot );
            }
            INSERT_TASK_zgetrf_panel_blocked_batched_flush( options, A, k, trsmUp,
                                                            Up(myrank), clargs, ipiv, pivot );

            assert( j <= minmn );

            /* Invalidate prevpiv because we don't need anymore */
            RUNTIME_pivot_invalidate( pivot, A->myrank, j-1 );

#if defined(CHAMELEON_USE_MPI)
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zpivot_allreduce( options, A, pivot, k, j, tempkn, ws->laswp );
#endif

            /*
             * Copy the row of the pivot in Up for the replicated TRSM.
             * If b is the last inner block then we don't have a TRSM.
             * If j+1 is a multiple of ib then we do the copy of the last
             * pivot row in the TRSM task.
             */
            if ( ( !last_ib ) && ( ( ( j + 1 ) % ib ) != 0 ) ) {
                INSERT_TASK_zgetrf_cpy_pivrow_in_Up( options, ws->Up, A->myrank,
                                                     k, j, ib, tempkn, pivot );
            }
        }
    }

    free( clargs );
    chameleon_data_flush( options->sequence, Up(myrank), request->flush );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, A->mt - k, ipiv, k );
    RUNTIME_pivot_flushone( options->sequence, pivot, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto( struct chameleon_pzgetrf_s *ws,
                               CHAM_desc_t                *A,
                               CHAM_ipiv_t                *ipiv,
                               CHAM_desc_pivot_t          *pivot,
                               int                         k,
                               RUNTIME_option_t           *options )
{
#if defined(CHAMELEON_USE_MPI)
    CHAM_reduce_t *reduce = &(ws->laswp->reduce);
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, k, reduce );
    if ( !reduce->involved ) {
        return;
    }
#endif

    /* TODO: Should be replaced by a function pointer */
    switch( ws->alg ) {
    case ChamGetrfPPivPerColumn:
        if ( ( ws->batch_size_blas2 > 0 ) || ws->batch_adaptive ) {
            chameleon_pzgetrf_panel_facto_percol_batched( ws, A, ipiv, pivot, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_percol( ws, A, ipiv, pivot, k, options );
        }
        break;

    case ChamGetrfPPiv:
    default:
        if ( ( ws->batch_size_blas2 > 0 ) || ws->batch_adaptive ) {
            chameleon_pzgetrf_panel_facto_blocked_batched( ws, A, ipiv, pivot, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_blocked( ws, A, ipiv, pivot, k, options );
        }
        break;
    }
}

#if defined(CHAMELEON_USE_MPI)
static inline void
chameleon_pzgetrf_panel_update_ws( struct chameleon_pzgetrf_s *ws,
                                   CHAM_desc_t                *A,
                                   int                         k,
                                   RUNTIME_option_t           *options )
{
    CHAM_context_t          *chamctxt = chameleon_context_self();
    const RUNTIME_request_t *request  = options->request;
    int m, n, tempkn, tempkm, p, involved, np;
    int lookahead = chamctxt->lookahead;
    int P         = chameleon_desc_datadist_get_iparam(A, 0);
    int Q         = chameleon_desc_datadist_get_iparam(A, 1);
    int lq        = (k % lookahead) * Q;
    int myp       = A->myrank / Q;

    for ( m = k+1; m < A->mt; m++ ) {
        if ( ( m % P ) != myp ) continue;

        chameleon_pzbcast_tile( ChamRowwise, ( k >= ws->ringswitch ) ? ChamBcastRing : ChamBcastFull,
                                A( m, k ), Wl( m, lq ), options );
        chameleon_data_flush( options->sequence, A(m, k), request->flush );
    }

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    np = P * Q;

    /* Send Akk for replicated trsm */
    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        for ( p = 0; p < np; p++ ) {
            involved = 0;
            for ( n = k+1; n < A->nt; n++ ) {
                if ( chameleon_p_involved_in_panelk_2dbc( A, n, p ) ) {
                    involved = 1;
                    break;
                }
            }
            if ( involved ) {
                INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempkn,
                                    A(k, k), Wu(p, k) );
            }
        }
    }
    else {
        involved = 0;
        for ( n = k+1; n < A->nt; n++ ) {
            if ( chameleon_involved_in_panelk_2dbc( A, n ) ) {
                involved = 1;
                break;
            }
        }
        if ( involved ) {
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempkn,
                                A(k, k), Wu(A->myrank, k) );
        }
    }
    chameleon_data_flush( options->sequence, A(k, k), request->flush );
}
#endif

void
chameleon_pzgetrf_panel_update( struct chameleon_pzgetrf_s *ws,
                                CHAM_desc_t                *A,
                                CHAM_ipiv_t                *ipiv,
                                int                         k,
                                int                         n,
                                RUNTIME_option_t           *options,
                                RUNTIME_sequence_t         *sequence )
{
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;
    CHAM_context_t             *chamctxt = chameleon_context_self();
    const RUNTIME_request_t    *request  = options->request;
    CHAM_reduce_t              *reduce   = &(ws->laswp->reduce);

    int m, tempkm, tempmm, tempnn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

    chameleon_pzlaswp_panel( ws->laswp, CHAMELEON_FALSE, ChamDirForward,
                             A, ipiv, k, n, options, sequence );

#if defined(CHAMELEON_USE_MPI)
    if ( RUNTIME_comm_size( chamctxt ) > 1 ) {
        int rankAmn;
        int lookahead = chamctxt->lookahead;
        int Q         = chameleon_desc_datadist_get_iparam(A, 1);
        int myq       = A->myrank % Q;
        int lq        = (k % lookahead) * Q;

        if ( reduce->involved ) {
            INSERT_TASK_ztrsm(
                options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, Wu(A->myrank, k),
                      Wu(A->myrank, n) );
        }

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            rankAmn = A->get_rankof( A, m, n );

            if ( A->myrank == rankAmn ) {
                INSERT_TASK_zgemm(
                    options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, Wl( m, myq + lq ),
                           Wu( A->myrank, n ),
                    zone,  A( m, n ) );
            }
        }

        if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                                Wu(A->myrank, n), A(k, n) );
        }
    }
    else
#endif
    {
        INSERT_TASK_ztrsm(
            options,
            ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
            tempkm, tempnn, A->mb,
            zone, A( k, k ),
                  Wu( A->myrank, n ) );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            INSERT_TASK_zgemm(
                options,
                ChamNoTrans, ChamNoTrans,
                tempmm, tempnn, A->mb, A->mb,
                mzone, A( m, k ),
                       Wu( A->myrank, n ),
                zone,  A( m, n ) );
        }

        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            Wu(A->myrank, n), A(k, n) );
    }

    chameleon_data_flush( options->sequence, Wu(A->myrank, n), request->flush );
    chameleon_data_flush( options->sequence, A(k, n), request->flush );
    (void)reduce;
    (void)chamctxt;
}

void
chameleon_pzgetrf_panel_permute_update( struct chameleon_pzgetrf_s *ws,
                                        CHAM_desc_t                *A,
                                        CHAM_ipiv_t                *ipiv,
                                        int                         k,
                                        int                         n,
                                        RUNTIME_option_t           *options,
                                        RUNTIME_sequence_t         *sequence )
{
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAM_context_t             *chamctxt = chameleon_context_self();
    const RUNTIME_request_t    *request  = options->request;
    CHAM_reduce_t              *reduce   = &(ws->laswp->reduce);
    void **clargs;
    int m, tempkm, tempmm, tempnn, withlacpy;

#if defined(CHAMELEON_USE_MPI)
    /* Initizalize the list of nodes invovlved in the panel n */
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, reduce );

    /* If on the rank who owns the ipiv array */
    if ( A->myrank == ipiv->get_rankof( ipiv, k, k ) ) {
        /* Exchange between all nodes involved the perm array */
        INSERT_TASK_zperm_allreduce_send_perm( options, ChamDirForward, ipiv, k, A->myrank,
                                               reduce->np_involved, reduce->proc_involved );

        /* Exchange between all nodes involved the invp array */
        INSERT_TASK_zperm_allreduce_send_invp_row( options, ChamDirForward, ipiv, k, A, k, n );
    }

    /* Bcast the top tile of the panel to all involved nodes */
    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank,
                                            reduce->np_involved, reduce->proc_involved );
    }

    /* If I'm not involved in the reduction, no need to go further */
    if ( !reduce->involved ) {
        return;
    }
#endif

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

    clargs = malloc( sizeof(char *) );
    *clargs = NULL;

    /* Copy tile A(k,n) into Wu(me,n) */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                        A(k, n), Wu(A->myrank, n) );
    options->withlacpy = withlacpy;

    /* Copy permuted rows of tiles A(m,n) into Wu(me,n) */
    for ( m = k; m < A->mt; m++ ) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        INSERT_TASK_zlaswp_get_batched( options, ChamLeft, ChamDirForward, m*A->mb, tempmm, tempnn, tempkm, (void *)ws->laswp, ipiv, k,
                                        A(m, n), Wu(A->myrank, n), clargs );
    }
    INSERT_TASK_zlaswp_get_batched_flush( options, ChamDirForward, ipiv, k, Wu(A->myrank, n), clargs );

#if defined(CHAMELEON_USE_MPI)
    /* Allreduce of Wu(me,n) to fully form the permuted tile A(k,n) and to replicate trsm */
    if ( ws->laswp->allreduce ) {
        INSERT_TASK_zperm_allreduce( options, ChamDirForward, A, Wu(A->myrank, n), ipiv, k, k, n, ws->laswp );
    }
    else {
        INSERT_TASK_zperm_reduce( options, ChamDirForward, A(k, n), ipiv, k, Wu(A->myrank, n), ws->laswp, A->myrank, n );
    }
#endif

#if defined(CHAMELEON_USE_MPI)
    int lookahead = chamctxt->lookahead;
    int Q         = chameleon_desc_datadist_get_iparam(A, 1);
    int myq       = A->myrank % Q;
    int lq        = (k % lookahead) * Q;

    /* Convert the tile ws(me,n) used to perform the allreduce into Wu(me,n) if needed */
    if ( reduce->np_involved != 1 ) {
        if ( reduce->alg_allreduce == ChamStarPUTasks ) {
            INSERT_TASK_zlaswp_ret( options, Ws(A->myrank, n), Wu(A->myrank, n) );
        }
        RUNTIME_perm_flush( sequence, A->myrank, Ws(A->myrank, n) );
    }

    if ( RUNTIME_comm_size( chamctxt ) > 1 ) {
        /* Trsm replicated on all processes involved in the column */
        if ( reduce->involved ) {
            INSERT_TASK_ztrsm(
                options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, Wu(A->myrank, k),
                      Wu(A->myrank, n) );
        }

        /* Copy the inverse permuted rows of the tile A(k,n) into A(m,n) and perform the gemm in the same task */
        for ( m = k + 1; m < A->mt; m++ ) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zlaswp_gemm( options, ChamLeft, ChamDirForward,
                                     m*A->mb, tempmm, tempnn, tempkm,
                                     (void *)ws->laswp, ipiv, k,
                                     A(k, n), Wl(m, myq + lq), Wu( A->myrank, n),
                                     A(m, n), clargs );

        }
        INSERT_TASK_zlaswp_gemm_flush( options, ChamDirForward, ipiv, k, A(k, n),
                                       Wu( A->myrank, n), clargs );

        /* Copy the tile Wu(me,n) into A(k,n) */
        if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                                Wu(A->myrank, n), A(k, n) );
        }
    }
    else
#endif
    {
        INSERT_TASK_ztrsm(
            options,
            ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
            tempkm, tempnn, A->mb,
            zone, A( k, k ),
                  Wu( A->myrank, n ) );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zlaswp_gemm( options, ChamLeft, ChamDirForward,
                                     m*A->mb, tempmm, tempnn, tempkm,
                                     (void *)ws->laswp, ipiv, k,
                                     A(k, n), A(m, k), Wu( A->myrank, n),
                                     A(m, n), clargs );
        }
        INSERT_TASK_zlaswp_gemm_flush( options, ChamDirForward, ipiv, k, A(k, n),
                                       Wu( A->myrank, n), clargs );

        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            Wu(A->myrank, n), A(k, n) );
    }

    free( clargs );
    chameleon_data_flush( options->sequence, Wu(A->myrank, n), request->flush );
    chameleon_data_flush( options->sequence, A(k, n), request->flush );
    (void)reduce;
    (void)chamctxt;
}

/**
 *  Parallel tile LU factorization with no pivoting - dynamic scheduling
 */
void chameleon_pzgetrf( struct chameleon_pzgetrf_s *ws,
                        CHAM_desc_t                *A,
                        CHAM_ipiv_t                *IPIV,
                        RUNTIME_sequence_t         *sequence,
                        RUNTIME_request_t          *request )
{
    CHAM_context_t    *chamctxt;
    RUNTIME_option_t   options;
    CHAM_desc_pivot_t *pivot = &(ws->pivot);

    int k, m, n;
    int min_mnt = chameleon_min( A->mt, A->nt );
    int kmin, kmax;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    kmin = chameleon_max( 0,       chamctxt->first_step );
    kmax = chameleon_min( min_mnt, chamctxt->last_step  );
    for (k = kmin; k < kmax; k++ ) {
        RUNTIME_iteration_push( chamctxt, k );

        options.priority = A->nt;
        /*
         * Do the panel factorization only if the current proc contributes in the
         * block column k.
         */
        options.forcesub = chameleon_involved_in_panelk_2dbc( A, k );
        if ( chameleon_involved_in_panelk_2dbc( A, k ) ) {
            chameleon_pzgetrf_panel_facto( ws, A, IPIV, pivot, k, &options );
        }
        options.forcesub = 0;

#if defined(CHAMELEON_USE_MPI)
        if ( RUNTIME_comm_size( chamctxt ) > 1 ) {
            chameleon_pzgetrf_panel_update_ws( ws, A, k, &options );
        }
#endif

        for (n = k+1; n < A->nt; n++) {
            options.priority = A->nt-n;
            if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                 chameleon_involved_in_panelk_2dbc( A, n ) )
            {
                ws->panel_permute_update( ws, A, IPIV, k, n, &options, sequence );
            }
        }

        /* Flush panel k */
        for (m = k+1; m < A->mt; m++) {
            chameleon_data_flush( sequence, A(m, k), request->flush );
        }
        chameleon_data_flush( sequence, Wu(A->myrank, k), request->flush );

        RUNTIME_iteration_pop( chamctxt );
    }
    if ( ws->Wl ) {
        CHAMELEON_Desc_Flush( ws->Wl, sequence );
    }
    CHAMELEON_Ipiv_Flush( IPIV, sequence );
    chameleon_pivot_destroy_submit( pivot, sequence );

    options.priority = 0;

    /* Backward pivoting */
    if ( ws->backperm_enabled ) {

        /* Disable allreduce for the backward permutation */
        ws->laswp->allreduce = 0;

        for (k = kmin+1; k < kmax; k++) {
            for (n = 0; n < k; n++) {
                if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                     chameleon_involved_in_panelk_2dbc( A, n ) )
                {
                    chameleon_pzlaswp_panel( ws->laswp, CHAMELEON_TRUE, ChamDirForward,
                                             A, IPIV, k, n, &options, sequence );
                }
                chameleon_data_flush( sequence, Wu(A->myrank, n), request->flush );
            }
            RUNTIME_ipiv_flushone( sequence, CHAMIPIV_PERM | CHAMIPIV_INVP, IPIV, k );
        }
    }
    CHAMELEON_Desc_Flush( ws->laswp->Wu, sequence );

    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize( &options, chamctxt );
}
