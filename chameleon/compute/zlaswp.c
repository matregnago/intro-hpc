/**
 *
 * @file zlaswp.c
 *
 * @copyright 2025-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp wrappers
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zlaswp_WS_Alloc - Allocate the required workspaces for laswp
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 *******************************************************************************
 *
 * @retval An allocated opaque pointer to use in CHAMELEON_laswp_Tile_Async()
 *         and to free with CHAMELEON_laswp_WS_Free().
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_zlaswp_WS_Free
 *
 */
void *
CHAMELEON_zlaswp_WS_Alloc( cham_side_t side, const CHAM_desc_t *A )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzlaswp_s *ws;
    CHAM_reduce_t              *reduce;
    int                         P            = chameleon_desc_datadist_get_iparam( A, 0 );
    int                         Q            = chameleon_desc_datadist_get_iparam( A, 1 );
    int                         max_involved = ( side == ChamLeft ) ? P : Q;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return NULL;
    }

    if ( A->get_rankof_init != chameleon_getrankof_2d ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_WS_Alloc", "zlaswp supports only 2DBC descriptors");
        return NULL;
    }
    ws = calloc( 1, sizeof(struct chameleon_pzlaswp_s) );

    ws->allreduce = 0;

    reduce = &(ws->reduce);
    reduce->alg_allreduce = ChamStarPUTasks;

#if defined (CHAMELEON_USE_MPI)
    reduce->proc_involved = malloc( sizeof( int ) * max_involved );
    reduce->involved      = 0;
    reduce->np_involved   = 0;
    reduce->arity         = chameleon_getenv_get_value_int( "CHAMELEON_ARITY", max_involved );
#else
    reduce->np_involved = 1;
#endif

    /*
     * Used only for testing purpose to switch from the reduce to allreduce algorithm in the tests
     */
    ws->allreduce = chameleon_getenv_get_value_int( "CHAMELEON_LASWP_ALLREDUCE", 0 );

    /*
     * Read the environment variable to define the allreduce mode to use
     */
    {
        char *allreduce = chameleon_getenv( "CHAMELEON_ALLREDUCE" );

        if ( allreduce != NULL ) {
            if ( strcasecmp( allreduce, "cham_spu_tasks" ) == 0 ) {
                reduce->alg_allreduce = ChamStarPUTasks;
            }
            else {
                chameleon_error( "CHAMELEON_zlaswp_WS_Alloc",
                                 "CHAMELEON_ALLREDUCE is not one of:\n"
                                 "   chameleon_starpu_tasks, chameleon_starpu,\n"
                                 "   chameleon_starpu_mpi, chameleon_mpi\n"
                                 "          => Switch back to chameleon_starpu_tasks\n" );
            }
        }
        chameleon_cleanenv( allreduce );
    }

    /*
     * Read the environment variable to setup the batch size
     */
    ws->batch_size_swap     = chameleon_getenv_get_value_int( "CHAMELEON_BATCH_SIZE", 0 );
    ws->batch_size_swap     = chameleon_getenv_get_value_int( "CHAMELEON_LASWP_BATCH_SIZE",     ws->batch_size_swap );
    ws->batch_size_swap_get = chameleon_getenv_get_value_int( "CHAMELEON_LASWP_BATCH_SIZE_GET", ws->batch_size_swap );
    ws->batch_size_swap_set = chameleon_getenv_get_value_int( "CHAMELEON_LASWP_BATCH_SIZE_SET", ws->batch_size_swap );
    if ( ws->batch_size_swap > CHAMELEON_BATCH_SIZE ) {
        chameleon_warning( "CHAMELEON_BATCH_SIZE",
                           "CHAMELEON_LASWP_BATCH_SIZE must be smaller than CHAMELEON_BATCH_SIZE:\n"
                           " please recompile with the right CHAMELEON_BATCH_SIZE, or reduce the CHAMELEON_LASWP_BATCH_SIZE value\n" );
        ws->batch_size_swap = CHAMELEON_BATCH_SIZE;
    }
    if ( ws->batch_size_swap_get > CHAMELEON_BATCH_SIZE ) {
        chameleon_warning( "CHAMELEON_BATCH_SIZE",
                           "CHAMELEON_LASWP_BATCH_SIZE_GET must be smaller than CHAMELEON_BATCH_SIZE:\n"
                           " please recompile with the right CHAMELEON_BATCH_SIZE, or reduce the CHAMELEON_LASWP_BATCH_SIZE_GET value\n" );
        ws->batch_size_swap_get = CHAMELEON_BATCH_SIZE;
    }
    if ( ws->batch_size_swap_set > CHAMELEON_BATCH_SIZE ) {
        chameleon_warning( "CHAMELEON_BATCH_SIZE",
                           "CHAMELEON_LASWP_BATCH_SIZE_SET must be smaller than CHAMELEON_BATCH_SIZE:\n"
                           " please recompile with the right CHAMELEON_BATCH_SIZE, or reduce the CHAMELEON_LASWP_BATCH_SIZE_SET value\n" );
        ws->batch_size_swap_set = CHAMELEON_BATCH_SIZE;
    }

    ws->Wu = malloc( sizeof(CHAM_desc_t) );
    if ( side == ChamLeft ) {
        chameleon_desc_init_2dtile( chamctxt, ws->Wu, "LASWP_Wu", ChamComplexDouble,
                                    A->mb, A->nb, A->mb * P * Q, A->n, P * Q, 1 );
    }
    else {
        chameleon_desc_init_2dtile( chamctxt, ws->Wu, "LASWP_Wu", ChamComplexDouble,
                                    A->mb, A->nb, A->m, A->nb * P * Q, 1, P * Q );
    }

    ws->ws.mt   = A->mt;
    ws->ws.nt   = A->nt;
    ws->ws.mb   = A->mb;
    ws->ws.nb   = A->nb;
    ws->ws.n    = A->n;
    ws->ws.m    = A->m;
    ws->ws.side = side;
    ws->ws.dtyp = A->dtyp;
    ws->ws.NP   = P * Q;

#if defined(CHAMELEON_USE_MPI)
    if ( reduce->alg_allreduce == ChamStarPUTasks ) {
        RUNTIME_perm_create( &(ws->ws) );
    }
#endif

    (void)max_involved;
    return ws;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Free the allocated workspaces for asynchronous laswp
 *
 *******************************************************************************
 *
 * @param[in,out] user_ws
 *          On entry, the opaque pointer allocated by
 *          CHAMELEON_zlaswp_WS_Alloc() On exit, all data are freed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_zlaswp_WS_Alloc
 *
 */
void
CHAMELEON_zlaswp_WS_Free( void *user_ws )
{
    struct chameleon_pzlaswp_s *ws = (struct chameleon_pzlaswp_s *)user_ws;

#if defined (CHAMELEON_USE_MPI)
    free( ws->reduce.proc_involved );
#endif

    if ( ws->Wu ) {
        chameleon_desc_destroy( ws->Wu );
        free( ws->Wu );
        ws->Wu = NULL;
    }
#if defined(CHAMELEON_USE_MPI)
    if ( ws->reduce.alg_allreduce == ChamStarPUTasks ) {
        RUNTIME_perm_destroy( &(ws->ws) );
    }
#endif

    free( ws );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Computes the permutation P*op(A) or op(A)*P where P is the permutation
 *         matrix generated from IPIV.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp_Tile
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_claswp
 * @sa CHAMELEON_dlaswp
 * @sa CHAMELEON_slaswp
 *
 */
int CHAMELEON_zlaswp( cham_side_t            side,
                      cham_dir_t             dir,
                      int                    M,
                      int                    N,
                      CHAMELEON_Complex64_t *A,
                      int                    LDA,
                      int                    K1,
                      int                    K2,
                      int                   *IPIV )
{
    int                 status;
    int                 NB;
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t         descAl, descAt;
    CHAM_ipiv_t         descIPIV;
    int                 K = ( side == ChamLeft ) ? M : N;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of M");
        return -2;
    }
    if ( N < 0 ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of N");
        return -3;
    }
    if ( LDA < chameleon_max( 1, M ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of LDA");
        return -5;
    }
    if ( ( K1 < 1 ) || ( K1 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return -6;
    }
    if ( ( K2 < 1 ) || ( K2 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return -7;
    }

    /* Quick return */
    if ( chameleon_min( N, M ) == 0 ) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zlaswp", "chameleon_tune() failed");
        return status;
    }

    /* Set NB */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    chameleon_ipiv_init( &descIPIV, side, descAt.mb, K2,
                         0, -1, -1,
                         1, 1, IPIV, NULL );

    /* Call the tile interface */
    ws = CHAMELEON_zlaswp_WS_Alloc( side, &descAt );
    CHAMELEON_zlaswp_Tile_Async( side, dir, &descAt, K1, K2, &descIPIV, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    CHAMELEON_zlaswp_WS_Free( ws );
    chameleon_ipiv_destroy( &descIPIV );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  @brief Tile equivalent of CHAMELEON_zlaswp().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_claswp_Tile
 * @sa CHAMELEON_dlaswp_Tile
 * @sa CHAMELEON_slaswp_Tile
 *
 */
int CHAMELEON_zlaswp_Tile( cham_side_t  side,
                           cham_dir_t   dir,
                           CHAM_desc_t *A,
                           int          K1,
                           int          K2,
                           CHAM_ipiv_t *IPIV )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    int                 K = ( side == ChamLeft ) ? A->m : A->n;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( ( K1 < 1 ) || ( K1 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( ( K2 < 1 ) || ( K2 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zlaswp_WS_Alloc( side, A );
    CHAMELEON_zlaswp_Tile_Async( side, dir, A, K1, K2, IPIV, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zlaswp_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  @brief Non-blocking equivalent of CHAMELEON_zlaswp_Tile().
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp
 * @sa CHAMELEON_zlaswp_Tile
 * @sa CHAMELEON_claswp_Tile_Async
 * @sa CHAMELEON_dlaswp_Tile_Async
 * @sa CHAMELEON_slaswp_Tile_Async
 *
 */
int CHAMELEON_zlaswp_Tile_Async( cham_side_t         side,
                                 cham_dir_t          dir,
                                 CHAM_desc_t        *A,
                                 int                 K1,
                                 int                 K2,
                                 CHAM_ipiv_t        *IPIV,
                                 void               *user_ws,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzlaswp_s *ws;
    RUNTIME_option_t            options;
    int                         k;
    int                         K = ( side == ChamLeft ) ? A->m : A->n;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( ( K1 < 1 ) || ( K1 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( ( K2 < 1 ) || ( K2 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if ( sequence->status == CHAMELEON_SUCCESS ) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if ( chameleon_desc_check(A) != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zlaswp_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( A->get_rankof_init != chameleon_getrankof_2d ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "zlaswp supports only 2DBC descriptors");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( A->mb != A->nb ) {
        chameleon_error("CHAMELEON_zlaswp_Tile_Async", "only matching tile sizes supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( chameleon_min( A->m, A->n ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( side, A );
    }
    else {
        ws = user_ws;
    }

    if ( IPIV->data != NULL ) {
        RUNTIME_options_init( &options, chamctxt, sequence, request );
        if ( side == ChamLeft ) {
            int tempkm, m0;

            for ( k = 0; k < IPIV->mt; k++ ) {
                tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                m0 = k * A->mb;
                if ( IPIV->get_rankof( IPIV, k, k ) != IPIV->myrank ) {
                    continue;
                }
                INSERT_TASK_ipiv_to_perm( &options, m0, tempkm, tempkm,
                                          K1 - 1, K2 - 1, A->mt - k, IPIV, k );
                RUNTIME_ipiv_flushone( sequence, CHAMIPIV_IPIV, IPIV, k );
            }
        }
        else {
            int tempkn, n0;

            for ( k = 0; k < IPIV->mt; k++ ) {
                tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                n0 = k * A->nb;
                if ( IPIV->get_rankof( IPIV, k, k ) != IPIV->myrank ) {
                    continue;
                }
                INSERT_TASK_ipiv_to_perm( &options, n0, tempkn, tempkn,
                                          K1 - 1, K2 - 1, A->nt - k, IPIV, k );
                RUNTIME_ipiv_flushone( sequence, CHAMIPIV_IPIV, IPIV, k );
            }
        }
    }

    if ( side == ChamLeft ) {
        chameleon_pzlaswp( ws, dir, A, IPIV, sequence, request );
    }
    else {
        chameleon_pzlaswpc( ws, dir, A, IPIV, sequence, request );
    }

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zlaswp_WS_Free( ws );
    }

    return CHAMELEON_SUCCESS;
}

