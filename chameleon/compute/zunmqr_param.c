/**
 *
 * @file zunmqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr_param wrappers
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
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zunmqr_param - Overwrites the general complex M-by-N matrix C with
 *
 *                  SIDE = 'L'     SIDE = 'R'
 *  TRANS = 'N':      Q * C          C * Q
 *  TRANS = 'C':      Q^H * C       C * Q^H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *        Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CHAMELEON_zgeqrf. Q is of order M if SIDE = ChamLeft
 *  and of order N if SIDE = ChamRight.
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
 *
 * @param[in] side
 *          Intended usage:
 *          = ChamLeft:  apply Q or Q^H from the left;
 *          = ChamRight: apply Q or Q^H from the right.
 *
 * @param[in] trans
 *          Intended usage:
 *          = ChamNoTrans:   no transpose, apply Q;
 *          = ChamConjTrans: conjugate transpose, apply Q^H.
 *
 * @param[in] M
 *          The number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix C. N >= 0.
 *
 * @param[in] K
 *          The number of elementary reflectors whose product defines
 *          the matrix Q.
 *          If side == ChamLeft,  M >= K >= 0.
 *          If side == ChamRight, N >= K >= 0.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by CHAMELEON_zgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          If side == ChamLeft,  LDA >= max(1,M).
 *          If side == ChamRight, LDA >= max(1,N).
 *
 * @param[in] descTS
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf.
 *
 * @param[in] descTT
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf.
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q^H*C.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zunmqr_param_Tile
 * @sa CHAMELEON_zunmqr_param_Tile_Async
 * @sa CHAMELEON_cunmqr
 * @sa CHAMELEON_dormqr
 * @sa CHAMELEON_sormqr
 * @sa CHAMELEON_zgeqrf
 *
 */
int CHAMELEON_zunmqr_param( const libhqr_tree_t *qrtree,
                        cham_side_t side, cham_trans_t trans, int M, int N, int K,
                        CHAMELEON_Complex64_t *A, int LDA,
                        CHAM_desc_t *descTS, CHAM_desc_t *descTT,
                        CHAMELEON_Complex64_t *C, int LDC )
{
    int NB, Am;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descCl, descCt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zunmqr_param", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ( side == ChamLeft ) {
        Am = M;
    } else {
        Am = N;
    }

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of side");
        return -1;
    }
    if ((trans != ChamConjTrans) && (trans != ChamNoTrans)){
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of trans");
        return -2;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of N");
        return -4;
    }
    if ((K < 0) || (K > Am)) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of K");
        return -5;
    }
    if (LDA < chameleon_max(1, Am)) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of LDA");
        return -7;
    }
    if (LDC < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zunmqr_param", "illegal value of LDC");
        return -10;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
    if (chameleon_min(M, chameleon_min(N, K)) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, K & N; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, K, N);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zunmqr_param", "chameleon_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
    NB   = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamLower,
                     A, NB, NB, LDA, K, Am, K, sequence, &request );
    chameleon_zlap2tile( chamctxt, "C", &descCl, &descCt, ChamDescInout, ChamUpperLower,
                     C, NB, NB, LDC, N, M,  N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zunmqr_param_Tile_Async( qrtree, side, trans, &descAt, descTS, descTT, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descTS, sequence );
    CHAMELEON_Desc_Flush( descTT, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descCl, &descCt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zunmqr_param_Tile - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by CHAMELEON_zgeqrf_Tile Q is of order M.
 *  All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = ChamLeft:  apply Q or Q^H from the left;
 *          = ChamRight: apply Q or Q^H from the right.
 *          Currently only ChamLeft is supported.
 *
 * @param[in] trans
 *          Intended usage:
 *          = ChamNoTrans:   no transpose, apply Q;
 *          = ChamConjTrans: conjugate transpose, apply Q^H.
 *          Currently only ChamConjTrans is supported.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by CHAMELEON_zgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by CHAMELEON_zgeqrf.
 *          Can be obtained with CHAMELEON_Alloc_Workspace_zgeqrf
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q^H*C.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zunmqr_param
 * @sa CHAMELEON_zunmqr_param_Tile_Async
 * @sa CHAMELEON_cunmqr_Tile
 * @sa CHAMELEON_dormqr_Tile
 * @sa CHAMELEON_sormqr_Tile
 * @sa CHAMELEON_zgeqrf_Tile
 *
 */
int CHAMELEON_zunmqr_param_Tile( const libhqr_tree_t *qrtree, cham_side_t side, cham_trans_t trans,
                             CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zunmqr_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zunmqr_param_Tile_Async( qrtree, side, trans, A, TS, TT, C, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( TS, sequence );
    CHAMELEON_Desc_Flush( TT, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  Non-blocking equivalent of CHAMELEON_zunmqr_param_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
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
 * @sa CHAMELEON_zunmqr_param
 * @sa CHAMELEON_zunmqr_param_Tile
 * @sa CHAMELEON_cunmqr_Tile_Async
 * @sa CHAMELEON_dormqr_Tile_Async
 * @sa CHAMELEON_sormqr_Tile_Async
 * @sa CHAMELEON_zgeqrf_Tile_Async
 *
 */
int CHAMELEON_zunmqr_param_Tile_Async( const libhqr_tree_t *qrtree,
                                   cham_side_t side, cham_trans_t trans,
                                   CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *C,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t D, *Dptr = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zunmqr_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zunmqr_param_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zunmqr_param_Tile", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zunmqr_param_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TS) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zunmqr_param_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TT) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zunmqr_param_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zunmqr_param_Tile", "invalid fourth descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || C->nb != C->mb) {
        chameleon_error("CHAMELEON_zunmqr_param_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((side != ChamLeft) && (side != ChamRight)) {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((trans != ChamConjTrans) && (trans != ChamNoTrans)){
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
    /*
     if (chameleon_min(M, chameleon_min(N, K)) == 0)
     return CHAMELEON_SUCCESS;
     */

#if defined(CHAMELEON_COPY_DIAG)
    {
        int n = chameleon_min( A->m, A->n );
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, A->m, n );
        Dptr = &D;
    }
#endif

    chameleon_pzunmqr_param( 1, qrtree, side, trans, A, C, TS, TT, Dptr, sequence, request );

    if ( Dptr != NULL ) {
        /* Flush data to be able to sync with runtimes such as PaRSEC */
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( C, sequence );
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
