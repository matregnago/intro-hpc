/**
 *
 * @file zgelqs_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqs_param wrappers
 *
 * @version 1.4.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgelqs_param - Compute a minimum-norm solution min || A*X - B || using the LQ factorization
 *  A = L*Q computed by CHAMELEON_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= M >= 0.
 *
 * @param[in] NRHS
 *          The number of columns of B. NRHS >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by CHAMELEON_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] descTS
 *          Auxiliary factorization data, computed by CHAMELEON_zgelqf.
 *
 * @param[in] descTT
 *          Auxiliary factorization data, computed by CHAMELEON_zgelqf.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= N.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgelqs_param_Tile
 * @sa CHAMELEON_zgelqs_param_Tile_Async
 * @sa CHAMELEON_cgelqs
 * @sa CHAMELEON_dgelqs
 * @sa CHAMELEON_sgelqs
 * @sa CHAMELEON_zgelqf
 *
 */
int CHAMELEON_zgelqs_param( const libhqr_tree_t *qrtree, int M, int N, int NRHS,
                        CHAMELEON_Complex64_t *A, int LDA,
                        CHAM_desc_t *descTS, CHAM_desc_t *descTT,
                        CHAMELEON_Complex64_t *B, int LDB )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqs_param", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgelqs_param", "illegal value of M");
        return -1;
    }
    if (N < 0 || M > N) {
        chameleon_error("CHAMELEON_zgelqs_param", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        chameleon_error("CHAMELEON_zgelqs_param", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgelqs_param", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, chameleon_max(1, N))) {
        chameleon_error("CHAMELEON_zgelqs_param", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqs_param", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInput, ChamUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, NRHS, N, NRHS, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgelqs_param_Tile_Async( qrtree, &descAt, descTS, descTT, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInput, ChamUpperLower, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descTS, sequence );
    CHAMELEON_Desc_Flush( descTT, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descBl, &descBt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgelqs_param_Tile - Computes a minimum-norm solution using previously computed
 *  LQ factorization.
 *  Tile equivalent of CHAMELEON_zgelqs_param().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by CHAMELEON_zgelqf.
 *
 * @param[in] TS
 *          Auxiliary factorization data, computed by CHAMELEON_zgelqf.
 *
 * @param[in] TT
 *          Auxiliary factorization data, computed by CHAMELEON_zgelqf.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgelqs_param
 * @sa CHAMELEON_zgelqs_param_Tile_Async
 * @sa CHAMELEON_cgelqs_Tile
 * @sa CHAMELEON_dgelqs_Tile
 * @sa CHAMELEON_sgelqs_Tile
 * @sa CHAMELEON_zgelqf_Tile
 *
 */
int CHAMELEON_zgelqs_param_Tile( const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                             CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqs_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgelqs_param_Tile_Async( qrtree, A, TS, TT, B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( TS, sequence );
    CHAMELEON_Desc_Flush( TT, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

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
 *  CHAMELEON_zgelqs_param_Tile_Async - Computes a minimum-norm solution using previously
 *  computed LQ factorization.
 *  Non-blocking equivalent of CHAMELEON_zgelqs_param_Tile().
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
 * @sa CHAMELEON_zgelqs_param
 * @sa CHAMELEON_zgelqs_param_Tile
 * @sa CHAMELEON_cgelqs_Tile_Async
 * @sa CHAMELEON_dgelqs_Tile_Async
 * @sa CHAMELEON_sgelqs_Tile_Async
 * @sa CHAMELEON_zgelqf_Tile_Async
 *
 */
int CHAMELEON_zgelqs_param_Tile_Async( const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_desc_t *subB;
    CHAM_desc_t *subA;
    CHAM_context_t *chamctxt;
    CHAM_desc_t D, *Dptr = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqs_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqs_param_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqs_param_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgelqs_param_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TS) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqs_param_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TT) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqs_param_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqs_param_Tile", "invalid fourth descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        chameleon_error("CHAMELEON_zgelqs_param_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
     return CHAMELEON_SUCCESS;
     }
     */
    /* subB = chameleon_desc_submatrix(B, A->m, 0, A->n-A->m, B->n);
     chameleon_pzlaset( ChamUpperLower, 0., 0., subB, sequence, request );
     free(subB); */

    subB = chameleon_desc_submatrix(B, 0, 0, A->m, B->n);
    subA = chameleon_desc_submatrix(A, 0, 0, A->m, A->m);
    chameleon_pztrsm( ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1.0, subA, subB, sequence, request );
    free(subA);
    free(subB);

#if defined(CHAMELEON_COPY_DIAG)
    {
        int m = chameleon_min( A->m, A->n );
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, m, A->n );
        Dptr = &D;
    }
#endif

    chameleon_pzunmlq_param( 1, qrtree, ChamLeft, ChamConjTrans, A, B, TS, TT, Dptr, sequence, request );
    if (Dptr != NULL) {
        /* Flush data to be able to sync with runtimes such as PaRSEC */
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Desc_Flush( B, sequence );
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
