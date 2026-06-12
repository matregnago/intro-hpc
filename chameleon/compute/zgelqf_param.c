/**
 *
 * @file zgelqf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf_param wrappers
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
 *  CHAMELEON_zgelqf_param - Computes the tile LQ factorization of a complex M-by-N matrix A: A = L * Q.
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
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data, required by CHAMELEON_zgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgelqf_param_Tile
 * @sa CHAMELEON_zgelqf_param_Tile_Async
 * @sa CHAMELEON_cgelqf
 * @sa CHAMELEON_dgelqf
 * @sa CHAMELEON_sgelqf
 * @sa CHAMELEON_zgelqs
 *
 */
int CHAMELEON_zgelqf_param( const libhqr_tree_t *qrtree, int M, int N,
                        CHAMELEON_Complex64_t *A, int LDA,
                        CHAM_desc_t *descTS, CHAM_desc_t *descTT )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqf_param", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zgelqf_param", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgelqf_param", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgelqf_param", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqf_param", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgelqf_param_Tile_Async( qrtree, &descAt, descTS, descTT, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descTS, sequence );
    CHAMELEON_Desc_Flush( descTT, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgelqf_param_Tile - Computes the tile LQ factorization of a matrix.
 *  Tile equivalent of CHAMELEON_zgelqf_param().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and below the diagonal of the array contain the m-by-min(M,N)
 *          lower trapezoidal matrix L (L is lower triangular if M <= N); the elements above the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by CHAMELEON_zgelqs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgelqf_param
 * @sa CHAMELEON_zgelqf_param_Tile_Async
 * @sa CHAMELEON_cgelqf_Tile
 * @sa CHAMELEON_dgelqf_Tile
 * @sa CHAMELEON_sgelqf_Tile
 * @sa CHAMELEON_zgelqs_Tile
 *
 */
int CHAMELEON_zgelqf_param_Tile( const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqf_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgelqf_param_Tile_Async( qrtree, A, TS, TT, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( TS, sequence );
    CHAMELEON_Desc_Flush( TT, sequence );

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
 *  CHAMELEON_zgelqf_param_Tile_Async - Computes the tile LQ factorization of a matrix.
 *  Non-blocking equivalent of CHAMELEON_zgelqf_param_Tile().
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
 * @sa CHAMELEON_zgelqf_param
 * @sa CHAMELEON_zgelqf_param_Tile
 * @sa CHAMELEON_cgelqf_Tile_Async
 * @sa CHAMELEON_dgelqf_Tile_Async
 * @sa CHAMELEON_sgelqf_Tile_Async
 * @sa CHAMELEON_zgelqs_Tile_Async
 *
 */
int CHAMELEON_zgelqf_param_Tile_Async( const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t D, *Dptr = NULL;
    int KT;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqf_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqf_param_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgelqf_param_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zgelqf_param_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TS) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqf_param_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TT) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgelqf_param_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zgelqf_param_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    /*
     if (chameleon_min(M, N) == 0)
     return CHAMELEON_SUCCESS;
     */

    if ( A->m < A->n ) {
        KT = A->mt;
    }
    else {
        KT = A->nt;
    }

#if defined(CHAMELEON_COPY_DIAG)
    {
        int m = chameleon_min( A->m, A->n );
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, m, A->n );
        Dptr = &D;
    }
#endif

    chameleon_pzgelqf_param( 1, KT, qrtree, A, TS, TT, Dptr, sequence, request );
    if (Dptr != NULL) {
        /* Flush data to be able to sync with runtimes such as PaRSEC */
        CHAMELEON_Desc_Flush( A, sequence );
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
