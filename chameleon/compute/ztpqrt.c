/**
 *
 * @file ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_ztpqrt - Computes a blocked QR factorization of a
 *  "triangular-pentagonal" matrix C, which is composed of a triangular block A
 *  and a pentagonal block B, using the compact representation for Q.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix B. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix B, and the order of the matrix
 *          A. N >= 0.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of B.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in,out] A
 *          On entry, the upper triangular N-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array
 *          contain the upper triangular matrix R.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
 *          are rectangular, and the last L rows are upper trapezoidal.
 *          On exit, B contains the pentagonal matrix V.  See Further Details.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.  LDB >= max(1,M).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data, required by CHAMELEON_zgeqrs to
 *          solve the system of equations, or by any function to apply the Q.
 *
 * @par Further Details:
 * =====================
 *
 *  The input matrix C is a (N+M)-by-N matrix
 *
 *               C = [ A ]
 *                   [ B ]
 *
 *  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
 *  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
 *  upper trapezoidal matrix B2:
 *
 *               B = [ B1 ]  <- (M-L)-by-N rectangular
 *                   [ B2 ]  <-     L-by-N upper trapezoidal.
 *
 *  The upper trapezoidal matrix B2 consists of the first L rows of a
 *  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
 *  B is rectangular M-by-N; if M=L=N, B is upper triangular.
 *
 *  The matrix W stores the elementary reflectors H(i) in the i-th column
 *  below the diagonal (of A) in the (N+M)-by-N input matrix C
 *
 *               C = [ A ]  <- upper triangular N-by-N
 *                   [ B ]  <- M-by-N pentagonal
 *
 *  so that W can be represented as
 *
 *               W = [ I ]  <- identity, N-by-N
 *                   [ V ]  <- M-by-N, same form as B.
 *
 *  Thus, all of information needed for W is contained on exit in B, which
 *  we call V above.  Note that V has the same form as B; that is,
 *
 *               V = [ V1 ] <- (M-L)-by-N rectangular
 *                   [ V2 ] <-     L-by-N upper trapezoidal.
 *
 *  The columns of V represent the vectors which define the H(i)'s.
 *
 *  The number of blocks is B = ceiling(N/NB), where each
 *  block is of order NB except for the last block, which is of order
 *  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
 *  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
 *  for the last block) T's are stored in the NB-by-N matrix T as
 *
 *               T = [T1 T2 ... TB].
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztpqrt_Tile
 * @sa CHAMELEON_ztpqrt_Tile_Async
 * @sa CHAMELEON_ctpqrt
 * @sa CHAMELEON_dtpqrt
 * @sa CHAMELEON_stpqrt
 * @sa CHAMELEON_zgeqrs
 *
 */
int CHAMELEON_ztpqrt( int M, int N, int L,
                  CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t *B, int LDB,
                  CHAM_desc_t *descT )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;
    CHAM_desc_t descBl, descBt;
    int minMN = chameleon_min( M, N );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztpqrt", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_ztpqrt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_ztpqrt", "illegal value of N");
        return -2;
    }
    if ((L < 0) || ((L > minMN) && (minMN > 0))) {
        chameleon_error("CHAMELEON_ztpqrt", "illegal value of N");
        return -3;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_ztpqrt", "illegal value of LDA");
        return -5;
    }
    if (LDB < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_ztpqrt", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (minMN == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune(CHAMELEON_FUNC_ZGELS, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztpqrt", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpper,
                     A, NB, NB, LDA, N, N, N, sequence, &request );
    chameleon_zlap2tile( chamctxt, "B", &descBl, &descBt, ChamDescInout, ChamUpperLower,
                     B, NB, NB, LDB, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_ztpqrt_Tile_Async( L, &descAt, &descBt, descT, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescInout, ChamUpper, sequence, &request );
    chameleon_ztile2lap( chamctxt, &descBl, &descBt,
                     ChamDescInout, ChamUpperLower, sequence, &request );
    CHAMELEON_Desc_Flush( descT, sequence );

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
 *  CHAMELEON_ztpqrt_Tile - Computes the tile QR factorization of a matrix.
 *  Tile equivalent of CHAMELEON_ztpqrt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by CHAMELEON_zgeqrs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_ztpqrt
 * @sa CHAMELEON_ztpqrt_Tile_Async
 * @sa CHAMELEON_ctpqrt_Tile
 * @sa CHAMELEON_dtpqrt_Tile
 * @sa CHAMELEON_stpqrt_Tile
 * @sa CHAMELEON_zgeqrs_Tile
 *
 */
int CHAMELEON_ztpqrt_Tile( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztpqrt_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_ztpqrt_Tile_Async( L, A, B, T, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( T, sequence );

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
 *  CHAMELEON_ztpqrt_Tile_Async - Computes the tile QR factorization of a matrix.
 *  Non-blocking equivalent of CHAMELEON_ztpqrt_Tile().
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
 * @sa CHAMELEON_ztpqrt
 * @sa CHAMELEON_ztpqrt_Tile
 * @sa CHAMELEON_ctpqrt_Tile_Async
 * @sa CHAMELEON_dtpqrt_Tile_Async
 * @sa CHAMELEON_stpqrt_Tile_Async
 * @sa CHAMELEON_zgeqrs_Tile_Async
 *
 */
int CHAMELEON_ztpqrt_Tile_Async( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_ztpqrt_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_ztpqrt_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_ztpqrt_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_ztpqrt_Tile", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(B) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztpqrt_Tile", "invalid second descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (chameleon_desc_check(T) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_ztpqrt_Tile", "invalid third descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_ztpqrt_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ( (L != 0) && (((B->m - L) % B->mb) != 0) ) {
        chameleon_error("CHAMELEON_ztpqrt_Tile", "Triangular part must be aligned with tiles");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* if (chamctxt->householder == ChamFlatHouseholder) { */
    chameleon_pztpqrt( L, A, B, T, sequence, request );
    /* } */
    /* else { */
    /*    chameleon_pztpqrtrh( A, T, CHAMELEON_RHBLK, sequence, request ); */
    /* } */

    return CHAMELEON_SUCCESS;
}
