/**
 *
 * @file zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe wrappers
 *
 * @version 1.4.0
 * @comment This file is a copy of zplghe.c
 *          wich has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Mathis Rade
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Generate a random hermitian (positive definite if 'bump' is large
 * enough) half-matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] uplo
 *          The part of the Matrix wich will be generated.
 *
 * @param[out] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa CHAMELEON_zplghe_Tile
 * @sa CHAMELEON_zplghe_Tile_Async
 * @sa CHAMELEON_cplghe
 * @sa CHAMELEON_dplghe
 * @sa CHAMELEON_splghe
 * @sa CHAMELEON_zplghe
 *
 */
int CHAMELEON_zplghe( double bump, cham_uplo_t uplo, int N,
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
        chameleon_fatal_error("CHAMELEON_zplghe", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        chameleon_error("CHAMELEON_zplghe", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zplghe", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_max(0, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zplghe", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescOutput, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zplghe_Tile_Async( bump, uplo, &descAt, seed, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescOutput, uplo, sequence, &request );

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
 * @brief Generate a random hermitian (positive definite if 'bump' is large
 * enough) half-matrix by tiles.
 *
 * Tile equivalent of CHAMELEON_zplghe().  Operates on matrices stored by
 * tiles.  All matrices are passed through descriptors.  All dimensions are
 * taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix.
 *
 * @param[in] uplo
 *          The part of the Matrix wich will be generated.
 *
 * @param[in] A
 *          On exit, The random hermitian matrix A generated.
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
 * @sa CHAMELEON_zplghe
 * @sa CHAMELEON_zplghe_Tile_Async
 * @sa CHAMELEON_cplghe_Tile
 * @sa CHAMELEON_dplghe_Tile
 * @sa CHAMELEON_splghe_Tile
 * @sa CHAMELEON_zplghe_Tile
 *
 */
int CHAMELEON_zplghe_Tile( double bump, cham_uplo_t uplo, CHAM_desc_t *A,
                           unsigned long long int seed )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zplghe_Tile_Async( bump, uplo, A, seed, sequence, &request );

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
 *  CHAMELEON_zplghe_Tile_Async - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 *  Non-blocking equivalent of CHAMELEON_zplghe_Tile().
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
 * @sa CHAMELEON_zplghe
 * @sa CHAMELEON_zplghe_Tile
 * @sa CHAMELEON_cplghe_Tile_Async
 * @sa CHAMELEON_dplghe_Tile_Async
 * @sa CHAMELEON_splghe_Tile_Async
 * @sa CHAMELEON_zplghe_Tile_Async
 *
 */
int CHAMELEON_zplghe_Tile_Async( double                 bump,
                                 cham_uplo_t            uplo,
                                 CHAM_desc_t           *A,
                                 unsigned long long int seed,
                                 RUNTIME_sequence_t    *sequence,
                                 RUNTIME_request_t     *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zplghe_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zplghe_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return CHAMELEON_SUCCESS;

    chameleon_pzplghe( bump, uplo, A, A->m, A->i, A->j, seed, sequence, request );

    return CHAMELEON_SUCCESS;
}
