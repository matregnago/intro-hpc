/**
 *
 * @file zplrnk.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnk wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
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
 *  CHAMELEON_zplrnk - Generate a random rank-k matrix by tiles.
 *
 *  Use the specific rank-k property of the following matrices multiplication
 *
 *  \f[ C = A \times B \f],
 *
 *  with A a random matrix of size m-by-k and B a random matrix of size k-by-n.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of C.
 *
 * @param[in] N
 *          The order of the matrix C. N >= 0.
 *
 * @param[in] K
 *          The rank of C.
 *
 * @param[out] C
 *           On exit, The random rank-k matrix C generated.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in] seedA
 *          The seed used in the random generation of matrix A.
 *
 * @param[in] seedB
 *          The seed used in the random generation of matrix B.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zplrnk_Tile
 * @sa CHAMELEON_zplrnk_Tile_Async
 * @sa CHAMELEON_cplrnk
 * @sa CHAMELEON_dplrnk
 * @sa CHAMELEON_splrnk
 *
 */
int CHAMELEON_zplrnk( int M, int N, int K,
                      CHAMELEON_Complex64_t *C, int LDC,
                      unsigned long long int seedA,
                      unsigned long long int seedB )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descCl, descCt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnk", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zplrnk", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zplrnk", "illegal value of N");
        return -2;
    }
    if (K < 0) {
        chameleon_error("CHAMELEON_zplrnk", "illegal value of K");
        return -3;
    }
    if (LDC < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zplrnk", "illegal value of LDC");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Set NT */
    NB = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "C", &descCl, &descCt, ChamDescOutput, ChamUpperLower,
                         C, NB, NB, LDC, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zplrnk_Tile_Async( K, &descCt, seedA, seedB, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descCl, &descCt,
                         ChamDescOutput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descCl, &descCt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zplrnk_Tile - Generate a random rank-k matrix by tiles.
 *  Tile equivalent of CHAMELEON_zplrnk().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] K
 *          The rank of C.
 *
 * @param[in] C
 *          On exit, The random rank-k matrix C generated.
 *
 * @param[in] seedA
 *          The seed used in the random generation of matrix A.
 *
 * @param[in] seedB
 *          The seed used in the random generation of matrix B.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zplrnk
 * @sa CHAMELEON_zplrnk_Tile_Async
 * @sa CHAMELEON_cplrnk_Tile
 * @sa CHAMELEON_dplrnk_Tile
 * @sa CHAMELEON_splrnk_Tile
 * @sa CHAMELEON_zplghe_Tile
 * @sa CHAMELEON_zplgsy_Tile
 *
 */
int CHAMELEON_zplrnk_Tile( int K, CHAM_desc_t *C,
                           unsigned long long int seedA,
                           unsigned long long int seedB)
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnk_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zplrnk_Tile_Async( K, C, seedA, seedB, sequence, &request );

    CHAMELEON_Desc_Flush( C, sequence );

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
 *  CHAMELEON_zplrnk_Tile_Async - Generate a random rank-k matrix by tiles.
 *  Non-blocking equivalent of CHAMELEON_zplrnk_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] K
 *          The rank of C.
 *
 * @param[in] C
 *          On exit, The random rank-k matrix C generated.
 *
 * @param[in] seedA
 *          The seed used in the random generation of matrix A.
 *
 * @param[in] seedB
 *          The seed used in the random generation of matrix B.
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
 * @sa CHAMELEON_zplrnk
 * @sa CHAMELEON_zplrnk_Tile
 * @sa CHAMELEON_cplrnk_Tile_Async
 * @sa CHAMELEON_dplrnk_Tile_Async
 * @sa CHAMELEON_splrnk_Tile_Async
 * @sa CHAMELEON_zplghe_Tile_Async
 * @sa CHAMELEON_zplgsy_Tile_Async
 *
 */
int CHAMELEON_zplrnk_Tile_Async( int K, CHAM_desc_t     *C,
                                 unsigned long long int seedA,
                                 unsigned long long int seedB,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnk_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnk_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnk_Tile", "NULL request");
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
    if (chameleon_desc_check(C) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zplrnk_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (C->nb != C->mb) {
        chameleon_error("CHAMELEON_zplrnk_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( C->m, C->n ) == 0)
        return CHAMELEON_SUCCESS;

    chameleon_pzplrnk( K, C, seedA, seedB, sequence, request );

    return CHAMELEON_SUCCESS;
}
