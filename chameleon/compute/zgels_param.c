/**
 *
 * @file zgels_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgels_param wrappers
 *
 * @version 1.4.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @author Alycia Lisito
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
 *  CHAMELEON_zgels_param - solves overdetermined or underdetermined linear
 *  systems involving an M-by-N matrix A using the QR or the LQ factorization of
 *  A.  It is assumed that A has full rank.  The following options are provided:
 *
 *  # trans = ChamNoTrans and M >= N: find the least squares solution of an
 *    overdetermined system, i.e., solve the least squares problem: minimize
 *    || B - A*X ||.
 *
 *  # trans = ChamNoTrans and M < N: find the minimum norm solution of an
 *    underdetermined system A * X = B.
 *
 *  Several right hand side vectors B and solution vectors X can be handled in a
 *  single call; they are stored as the columns of the M-by-NRHS right hand side
 *  matrix B and the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @param[in] qrtree
 *          The tree used for the factorization
 *
 * @param[in] trans
 *          Intended usage:
 *          = ChamNoTrans:   the linear system involves A;
 *          = ChamConjTrans: the linear system involves A^H.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the
 *          matrices B and X.  NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, if M >= N, A is overwritten by details of its QR
 *          factorization as returned by CHAMELEON_zgeqrf; if M < N, A is
 *          overwritten by details of its LQ factorization as returned by
 *          CHAMELEON_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descTS
 *          On exit, auxiliary factorization data.
 *
 * @param[out] descTT
 *          On exit, auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS matrix B of right hand side vectors, stored
 *          columnwise;
 *          On exit, if return value = 0, B is overwritten by the solution
 *          vectors, stored columnwise: if M >= N, rows 1 to N of B contain the
 *          least squares solution vectors; the residual sum of squares for the
 *          solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column; if M < N, rows 1 to N
 *          of B contain the minimum norm solution vectors;
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= MAX(1,M,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgels_param_Tile
 * @sa CHAMELEON_zgels_param_Tile_Async
 * @sa CHAMELEON_cgels
 * @sa CHAMELEON_dgels
 * @sa CHAMELEON_sgels
 *
 */
int CHAMELEON_zgels_param( const libhqr_tree_t *qrtree, cham_trans_t trans, int M, int N, int NRHS,
                           CHAMELEON_Complex64_t *A, int LDA,
                           CHAM_desc_t *descTS, CHAM_desc_t *descTT,
                           CHAMELEON_Complex64_t *B, int LDB )
{
    int i, j;
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgels_param", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (trans != ChamNoTrans) && (trans != ChamConjTrans) ) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of trans");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if (M < 0) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of N");
        return -3;
    }
    if (NRHS < 0) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of NRHS");
        return -4;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of LDA");
        return -6;
    }
    if (LDB < chameleon_max(1, chameleon_max(M, N))) {
        chameleon_error("CHAMELEON_zgels_param", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (chameleon_min(M, chameleon_min(N, NRHS)) == 0) {
        for (i = 0; i < chameleon_max(M, N); i++)
            for (j = 0; j < NRHS; j++)
                B[j*LDB+i] = 0.0;
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, N, NRHS);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgels_param", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    int maxMN = chameleon_max( M, N );

    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                         A, NB, NB, LDA, N,    M,     N,    sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                         B, NB, NB, LDB, NRHS, maxMN, NRHS, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgels_param_Tile_Async( qrtree, trans, &descAt, descTS, descTT, &descBt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, ChamUpperLower, sequence, &request );
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
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgels_param_Tile - solves overdetermined or underdetermined linear
 *  systems involving an M-by-N matrix A using the QR or the LQ factorization of
 *  A.
 *  Tile equivalent of CHAMELEON_zgels_param().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended usage:
 *          = ChamNoTrans:   the linear system involves A;
 *          = ChamConjTrans: the linear system involves A^H.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, if M >= N, A is overwritten by details of its QR
 *          factorization as returned by CHAMELEON_zgeqrf; if M < N, A is
 *          overwritten by details of its LQ factorization as returned by
 *          CHAMELEON_zgelqf.
 *
 * @param[out] TS
 *          On exit, auxiliary factorization data.
 *
 * @param[out] TT
 *          On exit, auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS matrix B of right hand side vectors, stored
 *          columnwise;
 *          On exit, if return value = 0, B is overwritten by the solution
 *          vectors, stored columnwise: if M >= N, rows 1 to N of B contain the
 *          least squares solution vectors; the residual sum of squares for the
 *          solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column; if M < N, rows 1 to N
 *          of B contain the minimum norm solution vectors;
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgels_param
 * @sa CHAMELEON_zgels_param_Tile_Async
 * @sa CHAMELEON_cgels_Tile
 * @sa CHAMELEON_dgels_Tile
 * @sa CHAMELEON_sgels_Tile
 *
 */
int CHAMELEON_zgels_param_Tile( const libhqr_tree_t *qrtree, cham_trans_t trans, CHAM_desc_t *A,
                                CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgels_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgels_param_Tile_Async( qrtree, trans, A, TS, TT, B, sequence, &request );

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
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zgels_param_Tile_Async - Solves overdetermined or underdetermined linear
 *  system of equations using the tile QR or the tile LQ factorization.
 *  Non-blocking equivalent of CHAMELEON_zgels_param_Tile().
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
 * @sa CHAMELEON_zgels_param
 * @sa CHAMELEON_zgels_param_Tile
 * @sa CHAMELEON_cgels_Tile_Async
 * @sa CHAMELEON_dgels_Tile_Async
 * @sa CHAMELEON_sgels_Tile_Async
 *
 */
int CHAMELEON_zgels_param_Tile_Async( const libhqr_tree_t *qrtree, cham_trans_t trans, CHAM_desc_t *A,
                                      CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *B,
                                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_desc_t *subA;
    CHAM_desc_t *subB;
    CHAM_context_t *chamctxt;
    CHAM_desc_t D, *Dptr = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgels_param_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgels_param_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgels_param_Tile", "NULL request");
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
    if ( (trans != ChamNoTrans) && (trans != ChamConjTrans) ) {
        chameleon_error("CHAMELEON_zgels_param_Tile_Async", "illegal value of trans");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgels_param_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TS) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgels_param_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(TT) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgels_param_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgels_param_Tile", "invalid fourth descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        chameleon_error("CHAMELEON_zgels_param_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    if (A->m >= A->n) {
#if defined(CHAMELEON_COPY_DIAG)
        {
            int n = chameleon_min( A->m, A->n );
            chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, A->m, n );
            Dptr = &D;
        }
#endif
        /*
         * compute QR factorization of A
         */
        chameleon_pzgeqrf_param( 1, A->nt, qrtree, A,
                                 TS, TT, Dptr, sequence, request );

        /* Perform the solve part */
        if ( trans == ChamNoTrans ) {
            /*
             * Least-Squares Problem min || A * X - B ||
             *
             * B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
             */
            chameleon_pzunmqr_param( 0, qrtree, ChamLeft, ChamConjTrans, A, B, TS, TT, Dptr, sequence, request );

            /*
             * B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
             */
            subB = chameleon_desc_submatrix(B, 0, 0, A->n, B->n);
            subA = chameleon_desc_submatrix(A, 0, 0, A->n, A->n);
            chameleon_pztrsm( ChamLeft, ChamUpper, ChamNoTrans, ChamNonUnit, 1.0, subA, subB, sequence, request );
        }
        else {
            /*
             * Underdetermined system of equations A**H * X = B
             *
             * B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)
             */
            subB = chameleon_desc_submatrix(B, 0, 0, A->n, B->n);
            subA = chameleon_desc_submatrix(A, 0, 0, A->n, A->n);
            chameleon_pztrsm( ChamLeft, ChamUpper, ChamConjTrans, ChamNonUnit, 1.0, subA, subB, sequence, request );

            /*
             * B(M+1:N,1:NRHS) = 0
             */
            /* TODO: enable subdescriptor with non aligned starting point in laset */
            /* free(subB); */
            /* subB = chameleon_desc_submatrix(B, A->n, 0, A->m-A->n, B->n); */
            /* chameleon_pzlaset( ChamUpperLower, 0., 0., subB, sequence, request ); */

            /*
             * B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)
             */
            chameleon_pzunmqr_param( 0, qrtree, ChamLeft, ChamNoTrans, A, B, TS, TT, Dptr, sequence, request );
        }
    }
    else {
#if defined(CHAMELEON_COPY_DIAG)
        {
            int m = chameleon_min( A->m, A->n );
            chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, m, A->n );
            Dptr = &D;
        }
#endif
        /*
         * Compute LQ factorization of A
         */
        chameleon_pzgelqf_param( 1, A->mt, qrtree, A, TS, TT, Dptr, sequence, request );

        /* Perform the solve part */
        if ( trans == ChamNoTrans ) {
            /*
             * Underdetermined system of equations A * X = B
             *
             * B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
             */
            subB = chameleon_desc_submatrix(B, 0, 0, A->m, B->n);
            subA = chameleon_desc_submatrix(A, 0, 0, A->m, A->m);
            chameleon_pztrsm( ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1.0, subA, subB, sequence, request );

            /*
             * B(M+1:N,1:NRHS) = 0
             */
            /* TODO: enable subdescriptor with non aligned starting point in laset */
            /* free(subB); */
            /* subB = chameleon_desc_submatrix(B, A->m, 0, A->n-A->m, B->n); */
            /* chameleon_pzlaset( ChamUpperLower, 0., 0., subB, sequence, request ); */

            /*
             * B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)
             */
            chameleon_pzunmlq_param( 0, qrtree, ChamLeft, ChamConjTrans, A, B, TS, TT, Dptr, sequence, request );
        }
        else {
            /*
             * Overdetermined system min || A**H * X - B ||
             *
             * B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
             */
            chameleon_pzunmlq_param( 0, qrtree, ChamLeft, ChamNoTrans, A, B, TS, TT, Dptr, sequence, request );

            /*
             * B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)
             */
            subB = chameleon_desc_submatrix(B, 0, 0, A->m, B->n);
            subA = chameleon_desc_submatrix(A, 0, 0, A->m, A->m);
            chameleon_pztrsm( ChamLeft, ChamLower, ChamConjTrans, ChamNonUnit, 1.0, subA, subB, sequence, request );
        }
    }

    free(subA);
    free(subB);

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
