/**
 *
 * @file zungqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr_param wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
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
 * @brief Generates an M-by-N matrix Q with orthonormal columns, which is defined
 * as the first N columns of a product of the elementary reflectors returned by
 * CHAMELEON_zgeqrf_param().
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
 *
 * @param[in] M
 *          The number of rows of the matrix Q. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix Q. N >= M.
 *
 * @param[in] K
 *          The number of columns of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by
 *          CHAMELEON_zgeqrf_param().
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] descTS
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf_param().
 *
 * @param[in] descTT
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf_param().
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zungqr_param_Tile
 * @sa CHAMELEON_zungqr_param_Tile_Async
 * @sa CHAMELEON_cungqr
 * @sa CHAMELEON_dorgqr
 * @sa CHAMELEON_sorgqr
 * @sa CHAMELEON_zgeqrf
 *
 */
int
CHAMELEON_zungqr_param( const libhqr_tree_t *qrtree, int M, int N, int K,
                        CHAMELEON_Complex64_t *A, int LDA,
                        CHAM_desc_t *descTS, CHAM_desc_t *descTT,
                        CHAMELEON_Complex64_t *Q, int LDQ )
{
    int                 NB;
    int                 status;
    CHAM_context_t *    chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t         descAl, descAt;
    CHAM_desc_t         descQl, descQt;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zungqr_param", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error( "CHAMELEON_zungqr_param", "illegal value of M" );
        return -1;
    }
    if ( ( N < 0 ) || ( N > M ) ) {
        chameleon_error( "CHAMELEON_zungqr_param", "illegal value of N" );
        return -2;
    }
    if ( ( K < 0 ) || ( K > N ) ) {
        chameleon_error( "CHAMELEON_zungqr_param", "illegal value of K" );
        return -3;
    }
    if ( LDA < chameleon_max( 1, M ) ) {
        chameleon_error( "CHAMELEON_zungqr_param", "illegal value of LDA" );
        return -5;
    }
    if ( LDQ < chameleon_max( 1, M ) ) {
        chameleon_error( "CHAMELEON_zungqr_param", "illegal value of LDQ" );
        return -8;
    }

    /* Quick return */
    if ( N <= 0 ) {
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = chameleon_tune( CHAMELEON_FUNC_ZGELS, M, N, 0 );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zungqr_param", "chameleon_tune() failed" );
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamLower,
                         A, NB, NB, LDA, N, M, K, sequence, &request );
    chameleon_zlap2tile( chamctxt, "Q", &descQl, &descQt, ChamDescInout, ChamUpperLower,
                         Q, NB, NB, LDQ, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zungqr_param_Tile_Async( qrtree, &descAt, descTS, descTT, &descQt,
                                       sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descQl, &descQt,
                         ChamDescInout, ChamUpperLower, sequence, &request );

    CHAMELEON_Desc_Flush( descTS, sequence );
    CHAMELEON_Desc_Flush( descTT, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descQl, &descQt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Generates an M-by-N matrix Q with orthonormal columns, which is defined
 * as the first N columns of a product of the elementary reflectors returned by
 * CHAMELEON_zgeqrf_param(). All matrices are passed through descriptors. All
 * dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by
 *          CHAMELEON_zgeqrf_param().
 *
 * @param[in] TS
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf_param().
 *
 * @param[in] TT
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf_param().
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zungqr_param
 * @sa CHAMELEON_zungqr_param_Tile_Async
 * @sa CHAMELEON_cungqr_Tile
 * @sa CHAMELEON_dorgqr_Tile
 * @sa CHAMELEON_sorgqr_Tile
 * @sa CHAMELEON_zgeqrf_Tile
 *
 */
int
CHAMELEON_zungqr_param_Tile( const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                             CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *Q )
{
    CHAM_context_t *    chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zungqr_param_Tile", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zungqr_param_Tile_Async( qrtree, A, TS, TT, Q, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( TS, sequence );
    CHAMELEON_Desc_Flush( TT, sequence );
    CHAMELEON_Desc_Flush( Q, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 * @brief Non-blocking equivalent of CHAMELEON_zungqr_param_Tile().
 *
 * This function may return before the computation is finished.
 * Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
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
 * @sa CHAMELEON_zungqr_param
 * @sa CHAMELEON_zungqr_param_Tile
 * @sa CHAMELEON_cungqr_Tile_Async
 * @sa CHAMELEON_dorgqr_Tile_Async
 * @sa CHAMELEON_sorgqr_Tile_Async
 * @sa CHAMELEON_zgeqrf_Tile_Async
 *
 */
int
CHAMELEON_zungqr_param_Tile_Async( const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                                   CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *Q,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t     D, *Dptr = NULL;
    int             M, N, K;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zungqr_param_Tile_Async", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zungqr_param_Tile_Async", "NULL sequence" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zungqr_param_Tile_Async", "NULL request" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if ( sequence->status == CHAMELEON_SUCCESS ) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED );
    }

    /* Check descriptors for correctness */
    if ( chameleon_desc_check( A ) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "invalid first descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( chameleon_desc_check( TS ) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "invalid second descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( chameleon_desc_check( TT ) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "invalid third descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( chameleon_desc_check( Q ) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "invalid fourth descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( ( A->nb != A->mb ) || ( Q->nb != Q->mb ) ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "only square tiles supported" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }

    M = Q->m;
    N = Q->n;
    K = chameleon_min( A->m, A->n );

    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "Incorrect value of M" );
        return chameleon_request_fail( sequence, request, -1 );
    }
    if ( N > M ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "Incorrect value of N" );
        return chameleon_request_fail( sequence, request, -2 );
    }
    if ( K > N ) {
        chameleon_error( "CHAMELEON_zungqr_param_Tile_Async", "Incorrect value of K" );
        return chameleon_request_fail( sequence, request, -3 );
    }

    /* Quick return */
    if ( N == 0 ) {
        return CHAMELEON_SUCCESS;
    }

#if defined( CHAMELEON_COPY_DIAG )
    {
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, A->m, K );
        Dptr = &D;
    }
#endif

    chameleon_pzlaset( ChamUpperLower, 0., 1., Q, sequence, request );
    chameleon_pzungqr_param( 1, chameleon_min( A->mt, A->nt ), qrtree,
                             A, Q, TS, TT, Dptr, sequence, request );

    if ( Dptr != NULL ) {
        /* Flush data to be able to sync with runtimes such as PaRSEC */
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( Q, sequence );
        CHAMELEON_Desc_Flush( TS, sequence );
        CHAMELEON_Desc_Flush( TT, sequence );
        CHAMELEON_Desc_Flush( Dptr, sequence );

        /* Wait and cleanup */
        chameleon_sequence_wait( chamctxt, sequence );
        chameleon_desc_destroy( Dptr );
    }
    (void)D;
    return CHAMELEON_SUCCESS;
}
