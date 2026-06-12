/**
 *
 * @file pzlaswpc.c
 *
 * @copyright 2025-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp parallel algorithm for column permutation.
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A,         m, n
#define Ws(m,n) &(ws->ws), m, n
#define Wu(m,n) ws->Wu,    m, n

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzlaswpc_panel_permute( struct chameleon_pzlaswp_s *ws,
                                  cham_dir_t                  dir,
                                  CHAM_desc_t                *A,
                                  CHAM_ipiv_t                *ipiv,
                                  int                         m,
                                  int                         k,
                                  RUNTIME_option_t           *options )
{
    int n;
    int tempmm, tempnn, tempkn;
    int withlacpy;

    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                        A(m, k), Wu(m, A->myrank) );
    options->withlacpy = withlacpy;

    INSERT_TASK_zlaswp_get( options, ChamRight, dir, k*A->nb, tempmm, tempkn, tempkn,
                            ipiv, k, A(m, k), Wu(m, A->myrank) );

    for ( n = k + 1; n < A->nt; n++ ) {
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        /* Extract selected rows into A(k, n) */
        INSERT_TASK_zlaswp_get( options, ChamRight, dir, n*A->nb, tempmm, tempnn, tempkn,
                                ipiv, k, A(m, n), Wu(m, A->myrank) );
        /* Copy rows from A(k,n) into their final position */
        INSERT_TASK_zlaswp_set( options, ChamRight, dir, n*A->nb, tempmm, tempnn, tempkn,
                                ipiv, k, A(m, k), A(m, n) );
    }

#if defined(CHAMELEON_USE_MPI)
    if ( ws->allreduce ) {
        INSERT_TASK_zperm_allreduce( options, dir, A, Wu(m, A->myrank), ipiv, k, m, k, ws );
    }
    else {
        INSERT_TASK_zperm_reduce( options, dir, A(m, k), ipiv, k, Wu(m, A->myrank), ws, m, A->myrank );
    }
#endif
}

/**
 *  Permutation of the panel n at step k
 *  This version batched sets of tasks together (one batch for the set, one for
 *  the set)
 */
static inline void
chameleon_pzlaswpc_panel_permute_batched( struct chameleon_pzlaswp_s *ws,
                                          cham_dir_t                  dir,
                                          CHAM_desc_t                *A,
                                          CHAM_ipiv_t                *ipiv,
                                          int                         m,
                                          int                         k,
                                          RUNTIME_option_t           *options )
{
    int n;
    int tempmm, tempnn, tempkn;
    int withlacpy;

    void **clargs = malloc( sizeof(char *) );
    *clargs = NULL;

    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                        A(m, k), Wu(m, A->myrank) );
    options->withlacpy = withlacpy;

    for ( n = k; n < A->nt; n++ ) {
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        INSERT_TASK_zlaswp_get_batched( options, ChamRight, dir, n*A->nb, tempmm, tempnn, tempkn, (void *)ws, ipiv, k,
                                        A(m, n), Wu(m, A->myrank), clargs );
    }
    INSERT_TASK_zlaswp_get_batched_flush( options, dir, ipiv, k, Wu(m, A->myrank), clargs );

    for ( n = k + 1; n < A->nt; n++ ) {
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        INSERT_TASK_zlaswp_set_batched( options, ChamRight, dir, n*A->nb, tempmm, tempnn, tempkn, (void *)ws, ipiv, k,
                                        A(m, k), A(m, n), clargs );
    }
    INSERT_TASK_zlaswp_set_batched_flush( options, dir, ipiv, k, A(m, k), clargs );

#if defined(CHAMELEON_USE_MPI)
    if ( ws->allreduce ) {
        INSERT_TASK_zperm_allreduce( options, dir, A, Wu(m, A->myrank), ipiv, k, m, k, ws );
    }
    else {
        INSERT_TASK_zperm_reduce( options, dir, A(m, k), ipiv, k, Wu(m, A->myrank), ws, m, A->myrank );
    }
#endif

    free( clargs );
}

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzlaswpc_panel_permute_batched2( struct chameleon_pzlaswp_s *ws,
                                           cham_dir_t                  dir,
                                           CHAM_desc_t                *A,
                                           CHAM_ipiv_t                *ipiv,
                                           int                         m,
                                           int                         k,
                                           RUNTIME_option_t           *options )
{
    int n;
    int tempmm, tempnn, tempkn;
    int withlacpy;

    void **clargs = malloc( sizeof(char *) );
    *clargs = NULL;

    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                        A(m, k), Wu(m, A->myrank) );
    options->withlacpy = withlacpy;

    INSERT_TASK_zlaswp_get( options, ChamRight, dir, k*A->nb, tempmm, tempkn, tempkn,
                            ipiv, k, A(m, k), Wu(m, A->myrank) );

    for ( n = k + 1; n < A->nt; n++ ) {
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        INSERT_TASK_zlaswp_batched( options, ChamRight, dir, n*A->nb, tempmm, tempnn, tempkn, (void *)ws, ipiv, k,
                                    A(m, n), A(m, k), Wu(m, A->myrank), clargs );
    }
    INSERT_TASK_zlaswp_batched_flush( options, dir, ipiv, k, A(m, k), Wu(m, A->myrank), clargs );

#if defined(CHAMELEON_USE_MPI)
    if ( ws->allreduce ) {
        INSERT_TASK_zperm_allreduce( options, dir, A, Wu(m, A->myrank), ipiv, k, m, k, ws );
    }
    else {
        INSERT_TASK_zperm_reduce( options, dir, A(m, k), ipiv, k, Wu(m, A->myrank), ws, m, A->myrank );
    }
#endif

    free( clargs );
}

static inline void
chameleon_pzlaswpc_panel( struct chameleon_pzlaswp_s *ws,
                          cham_bool_t                 inplace,
                          cham_dir_t                  dir,
                          CHAM_desc_t                *A,
                          CHAM_ipiv_t                *ipiv,
                          int                         m,
                          int                         k,
                          RUNTIME_option_t           *options,
                          RUNTIME_sequence_t         *sequence )
{
    const RUNTIME_request_t *request = options->request;
    CHAM_reduce_t           *reduce  = &(ws->reduce);
    int                      tempmm, tempkn;

#if defined(CHAMELEON_USE_MPI)
    /* Initizalize the list of nodes invovlved in the row panel m */
    chameleon_get_proc_involved_in_rowpanelk_2dbc( A, m, k, reduce );

    /* If on the rank who owns the ipiv array */
    if ( A->myrank == ipiv->get_rankof( ipiv, k, k ) ) {
        /* Exchange between all nodes involved the perm array */
        INSERT_TASK_zperm_allreduce_send_perm( options, dir, ipiv, k, A->myrank,
                                               reduce->np_involved, reduce->proc_involved );

        /* Exchange between all nodes involved the invp array */
        INSERT_TASK_zperm_allreduce_send_invp_col( options, dir, ipiv, k, A, m, k );
    }

    /* Bcast the top tile of the panel to all involved nodes */
    if ( A->myrank == chameleon_getrankof_2d( A, m, k ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, m, k, A->myrank,
                                            reduce->np_involved, reduce->proc_involved );
    }

    /* If I'm not involved in the reduction, no need to go further */
    if ( !reduce->involved ) {
        return;
    }
#endif

    /*
     * Perform the permutation on the panel, the final top tile is stored in
     * Wu(m,...) or Ws(m,...) when done depending on the configuration.
     */
    if ( ws->batch_size_swap == 0 ){
        chameleon_pzlaswpc_panel_permute( ws, dir, A, ipiv, m, k, options );
    }
    else {
        chameleon_pzlaswpc_panel_permute_batched( ws, dir, A, ipiv, m, k, options );
    }

    /*
     * Now, we need to set the final top tile to the right position. Either
     * directly in A (when peforming a standalone swap operation), or in the
     * temporary replicated buffer to perform a replicated trsm (LU
     * factorization for example).
     */
    if ( inplace ) {
        if ( ws->reduce.np_involved == 1 ) {
            /*
             * If we just perform a laswp, A must be equal to Atop, then we copy the
             * final version of the tile to its final location.
             */
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                                Wu(m, A->myrank), A( m, k ) );
        }
#if defined(CHAMELEON_USE_MPI)
        else {
            if ( reduce->alg_allreduce == ChamStarPUTasks ) {
                /*
                 * Let's copy the final version of A to its final position
                 */
                INSERT_TASK_zlaswp_ret( options, Ws(m, A->myrank), A(m, k) );
                RUNTIME_perm_flush( sequence, A->myrank, Ws(m, A->myrank) );
            }
        }
#endif
        chameleon_data_flush( sequence, A(m, k), request->flush );
    }
#if defined(CHAMELEON_USE_MPI)
    else { /* Outofplace with replication */
        if ( reduce->np_involved != 1 ) {
            if ( reduce->alg_allreduce == ChamStarPUTasks ) {
                INSERT_TASK_zlaswp_ret( options, Ws(m, A->myrank), Wu(m, A->myrank) );
            }
            RUNTIME_perm_flush( sequence, A->myrank, Ws(m, A->myrank) );
        }
    }
#endif
    (void)reduce;
}

void
chameleon_pzlaswpc( struct chameleon_pzlaswp_s *ws,
                    cham_dir_t                  dir,
                    CHAM_desc_t                *A,
                    CHAM_ipiv_t                *IPIV,
                    RUNTIME_sequence_t         *sequence,
                    RUNTIME_request_t          *request )
{
    CHAM_context_t   *chamctxt;
    RUNTIME_option_t  options;

    int m, k;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    if ( dir == ChamDirForward ) {
        for ( k = 0; k < IPIV->mt; k++ ) {
            for ( m = 0; m < A->mt; m++ ) {
                options.priority = A->mt-m;

                chameleon_pzlaswpc_panel( ws, CHAMELEON_TRUE, dir, A, IPIV, m, k, &options, sequence );
            }
            RUNTIME_ipiv_flushone( sequence, CHAMIPIV_PERM | CHAMIPIV_INVP, IPIV, k );
        }
    }
    else {
        for ( k = IPIV->mt - 1; k > -1; k-- ) {
            for ( m = 0; m < A->mt; m++ ) {
                options.priority = A->mt-m;
                chameleon_pzlaswpc_panel( ws, CHAMELEON_TRUE, dir, A, IPIV, m, k, &options, sequence );
            }
            RUNTIME_ipiv_flushone( sequence, CHAMIPIV_PERM | CHAMIPIV_INVP, IPIV, k );
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize( &options, chamctxt );
}
