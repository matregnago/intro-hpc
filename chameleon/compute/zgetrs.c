/**
 *
 * @file zgetrs.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrs wrappers
 *
 * @version 1.4.0
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
 *  @brief Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization with partial pivoting computed by CHAMELEON_zgetrf.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = ChamNoTrans:   A * X = B     (No transpose)
 *          = ChamTrans:     A^T * X = B  (Transpose)
 *          = ChamConjTrans: A^H * X = B  (Conjugate transpose)
 *          Only ChamNoTrans and ChamTrans are supported.
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by CHAMELEON_zgetrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] IPIV
 *          On entry, ipiv descriptor associated to A and created with
 *          CHAMELEON_Ipiv_Create().
 *          On exit, it contains the pivot indices associated to the PLU
 *          factorization of A.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrs_Tile
 * @sa CHAMELEON_zgetrs_Tile_Async
 * @sa CHAMELEON_cgetrs
 * @sa CHAMELEON_dgetrs
 * @sa CHAMELEON_sgetrs
 * @sa CHAMELEON_zgetrf
 *
 */
int CHAMELEON_zgetrs( cham_trans_t trans, int N, int NRHS,
                      CHAMELEON_Complex64_t *A, int LDA,
                      int *IPIV,
                      CHAMELEON_Complex64_t *B, int LDB )
{
    int                         NB;
    int                         status;
    CHAM_context_t             *chamctxt;
    CHAM_ipiv_t                 descIPIV;
    RUNTIME_sequence_t         *sequence = NULL;
    RUNTIME_request_t           request  = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t                 descAl, descAt;
    CHAM_desc_t                 descBl, descBt;
    struct chameleon_pzgetrf_s *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zgetrs", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( ( trans != ChamTrans ) && ( trans != ChamNoTrans ) ) {
        chameleon_error("CHAMELEON_zgetrs", "Only ChamTrans and ChamNoTrans are supported");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( N < 0 ) {
        chameleon_error("CHAMELEON_zgetrs", "illegal value of N");
        return -2;
    }
    if ( NRHS < 0 ) {
        chameleon_error("CHAMELEON_zgetrs", "illegal value of NRHS");
        return -3;
    }
    if ( LDA < chameleon_max( 1, N ) ) {
        chameleon_error("CHAMELEON_zgetrs", "illegal value of LDA");
        return -5;
    }
    if ( LDB < chameleon_max( 1, N ) ) {
        chameleon_error("CHAMELEON_zgetrs", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if ( chameleon_min( N, NRHS ) == 0 )
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = chameleon_tune( CHAMELEON_FUNC_ZGESV, N, N, NRHS );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zgetrs", "chameleon_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                         B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    ws = CHAMELEON_zlaswp_WS_Alloc( ChamLeft, &descBt );
    chameleon_ipiv_init( &descIPIV, ChamLeft, descAt.mb, N,
                         0, -1, -1,
                         1, 1, IPIV, NULL );

    /* Call the tile interface */
    CHAMELEON_zgetrs_Tile_Async( trans, &descAt, &descIPIV, &descBt, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                         ChamDescInout, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ipiv_destroy( &descIPIV );
    CHAMELEON_zlaswp_WS_Free( ws );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  @brief Solves a system of linear equations using previously
 *  computed LU factorization with partial pivoting.
 *  Tile equivalent of CHAMELEON_zgetrs().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = ChamNoTrans:   A * X = B     (No transpose)
 *          = ChamTrans:     A^T * X = B  (Transpose)
 *          = ChamConjTrans: A^H * X = B  (Conjugate transpose)
 *          Only ChamNoTrans and ChamTrans are supported.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by CHAMELEON_zgetrf.
 *
 * @param[in] IPIV
 *          On entry, ipiv descriptor associated to A and created with
 *          CHAMELEON_Ipiv_Create().
 *          On exit, it contains the pivot indices associated to the PLU
 *          factorization of A.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrs
 * @sa CHAMELEON_zgetrs_Tile_Async
 * @sa CHAMELEON_cgetrs_Tile
 * @sa CHAMELEON_dgetrs_Tile
 * @sa CHAMELEON_sgetrs_Tile
 * @sa CHAMELEON_zgetrf_Tile
 *
 */
int CHAMELEON_zgetrs_Tile( cham_trans_t trans,
                           CHAM_desc_t *A,
                           CHAM_ipiv_t *IPIV,
                           CHAM_desc_t *B )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( ( trans != ChamTrans ) && ( trans != ChamNoTrans ) ) {
        chameleon_error("CHAMELEON_zgetrs", "Only ChamTrans and ChamNoTrans are supported");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zgetrs_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zlaswp_WS_Alloc( ChamLeft, B );

    CHAMELEON_zgetrs_Tile_Async( trans, A, IPIV, B, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

    CHAMELEON_zlaswp_WS_Free( ws );

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
 *  @brief Solves a system of linear equations using previously
 *  computed LU factorization with partial pivoting.
 *  Non-blocking equivalent of CHAMELEON_zgetrs_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = ChamNoTrans:   A * X = B     (No transpose)
 *          = ChamTrans:     A^T * X = B  (Transpose)
 *          = ChamConjTrans: A^H * X = B  (Conjugate transpose)
 *          Only ChamNoTrans and ChamTrans are supported.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] IPIV
 *          On entry, ipiv descriptor associated to A and created with
 *          CHAMELEON_Ipiv_Create().
 *          On exit, it contains the pivot indices associated to the PLU
 *          factorization of A.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the N-by-NRHS solution matrix X.
 *
 * @param[in,out] user_ws
 *          The opaque pointer to pre-allocated getrf workspace through
 *          CHAMELEON_zgetrf_WS_Alloc() for B. If user_ws is NULL, it is automatically
 *          allocated, but BE CAREFULL as it switches the call from asynchronous
 *          to synchronous call.*
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
 * @sa CHAMELEON_zgetrs
 * @sa CHAMELEON_zgetrs_Tile
 * @sa CHAMELEON_cgetrs_Tile_Async
 * @sa CHAMELEON_dgetrs_Tile_Async
 * @sa CHAMELEON_sgetrs_Tile_Async
 * @sa CHAMELEON_zgetrf_Tile_Async
 *
 */
int CHAMELEON_zgetrs_Tile_Async( cham_trans_t        trans,
                                 CHAM_desc_t        *A,
                                 CHAM_ipiv_t        *IPIV,
                                 CHAM_desc_t        *B,
                                 void               *user_ws,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzlaswp_s *ws;
    RUNTIME_option_t            options;
    int                         k, tempkm;
    CHAMELEON_Complex64_t       zone = (CHAMELEON_Complex64_t)1.;

    chamctxt = chameleon_context_self();
    if ( ( trans != ChamTrans ) && ( trans != ChamNoTrans ) ) {
        chameleon_error("CHAMELEON_zgetrs", "Only ChamTrans and ChamNoTrans are supported");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zgetrs_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error("CHAMELEON_zgetrs_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error("CHAMELEON_zgetrs_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgetrs_Tile", "invalid first descriptor");
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( chameleon_desc_check( B ) != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zgetrs_Tile", "invalid third descriptor");
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    /* Check input arguments */
    if ( ( A->nb != A->mb ) || ( B->nb != B->mb ) ) {
        chameleon_error("CHAMELEON_zgetrs_Tile", "only square tiles supported");
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( ChamLeft, B );
    }
    else {
        ws = user_ws;
    }

    if ( IPIV->data != NULL ) {
        RUNTIME_options_init( &options, chamctxt, sequence, request );
        for ( k = 0; k < A->mt; k++ ) {
            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            INSERT_TASK_ipiv_to_perm( &options, k * A->mb, tempkm, tempkm,
                                      0, A->m, A->mt - k,
                                      IPIV, k );
        }
        chameleon_sequence_wait( chamctxt, sequence );
    }

    if ( trans == ChamNoTrans ) {
        chameleon_pzlaswp( ws, ChamDirForward,
                           B, IPIV, sequence, request );

        chameleon_pztrsm( ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                          zone, A, B, sequence, request );
        chameleon_pztrsm( ChamLeft, ChamUpper, ChamNoTrans, ChamNonUnit,
                          zone, A, B, sequence, request );
    }
    else {
        chameleon_pztrsm( ChamLeft, ChamUpper, ChamTrans, ChamNonUnit,
                          zone, A, B, sequence, request );
        chameleon_pztrsm( ChamLeft, ChamLower, ChamTrans, ChamUnit,
                          zone, A, B, sequence, request );

        chameleon_pzlaswp( ws, ChamDirBackward,
                           B, IPIV, sequence, request );
    }

    if ( user_ws == NULL ) {
        CHAMELEON_zlaswp_WS_Free( ws );
    }

    return CHAMELEON_SUCCESS;
}
