/**
 *
 * @file zplrnt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt wrappers
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
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
 *  CHAMELEON_zplrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zplrnt_Tile
 * @sa CHAMELEON_zplrnt_Tile_Async
 * @sa CHAMELEON_cplrnt
 * @sa CHAMELEON_dplrnt
 * @sa CHAMELEON_splrnt
 * @sa CHAMELEON_zplghe
 * @sa CHAMELEON_zplgsy
 *
 */
int CHAMELEON_zplrnt( int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA,
                  unsigned long long int seed )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnt", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zplrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zplrnt", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zplrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zplrnt", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescOutput, ChamUpperLower,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zplrnt_Tile_Async( &descAt, seed, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescOutput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zplrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of CHAMELEON_zplrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zplrnt
 * @sa CHAMELEON_zplrnt_Tile_Async
 * @sa CHAMELEON_cplrnt_Tile
 * @sa CHAMELEON_dplrnt_Tile
 * @sa CHAMELEON_splrnt_Tile
 * @sa CHAMELEON_zplghe_Tile
 * @sa CHAMELEON_zplgsy_Tile
 *
 */
int CHAMELEON_zplrnt_Tile( CHAM_desc_t *A,
                       unsigned long long int seed )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnt_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zplrnt_Tile_Async( A, seed, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

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
 *  CHAMELEON_zplrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of CHAMELEON_zplrnt_Tile().
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
 * @sa CHAMELEON_zplrnt
 * @sa CHAMELEON_zplrnt_Tile
 * @sa CHAMELEON_cplrnt_Tile_Async
 * @sa CHAMELEON_dplrnt_Tile_Async
 * @sa CHAMELEON_splrnt_Tile_Async
 * @sa CHAMELEON_zplghe_Tile_Async
 * @sa CHAMELEON_zplgsy_Tile_Async
 *
 */
int CHAMELEON_zplrnt_Tile_Async( CHAM_desc_t     *A,
                             unsigned long long int seed,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnt_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnt_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zplrnt_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zplrnt_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return CHAMELEON_SUCCESS;

    chameleon_pzplrnt( A, A->m, A->i, A->j, seed, sequence, request );

    return CHAMELEON_SUCCESS;
}
