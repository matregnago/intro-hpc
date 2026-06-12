/**
 *
 * @file zplgtr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgtr wrappers
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Generate a random trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix if needed.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper trapezoidal or lower trapezoidal:
 *          = ChamUpper: Only the upper trapezoidal part of A is generated;
 *          = ChamLower: Only the lower trapezoidal part of A is generated;
 *          = ChamUpperLower: the full A is generated. Note that in that case
 *          the matrix is symmetric, and it is thus not recommended to call
 *          this function.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, the random trapezoidal matrix A generated.
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
 * @sa CHAMELEON_zplgtr_Tile
 * @sa CHAMELEON_zplgtr_Tile_Async
 * @sa CHAMELEON_cplgtr
 * @sa CHAMELEON_dplgtr
 * @sa CHAMELEON_splgtr
 * @sa CHAMELEON_zplgtr
 *
 */
int CHAMELEON_zplgtr( CHAMELEON_Complex64_t bump, cham_uplo_t uplo,
                      int M, int N, CHAMELEON_Complex64_t *A, int LDA,
                      unsigned long long int seed )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zplgtr", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( N < 0 ) {
        chameleon_error("CHAMELEON_zplgtr", "illegal value of N");
        return -2;
    }
    if ( M < 0 ) {
        chameleon_error("CHAMELEON_zplgtr", "illegal value of M");
        return -3;
    }
    if ( LDA < chameleon_max(1, M) ) {
        chameleon_error("CHAMELEON_zplgtr", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if ( chameleon_max(0, N) == 0 || chameleon_max(0, M) == 0 )
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zplgtr", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescOutput, uplo,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zplgtr_Tile_Async( bump, uplo, &descAt, seed, sequence, &request );

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
 * @brief Generate a random trapezoidal matrix. Tile equivalent of CHAMELEON_zplgtr().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix if needed.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper trapezoidal or lower trapezoidal:
 *          = ChamUpper: Only the upper trapezoidal part of A is generated;
 *          = ChamLower: Only the lower trapezoidal part of A is generated;
 *          = ChamUpperLower: the full A is generated. Note that in that case
 *          the matrix is symmetric, and it is thus not recommended to call
 *          this function.
 *
 * @param[out] A
 *          On exit, the random trapezoidal matrix A generated.
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
 * @sa CHAMELEON_zplgtr
 * @sa CHAMELEON_zplgtr_Tile_Async
 * @sa CHAMELEON_cplgtr_Tile
 * @sa CHAMELEON_dplgtr_Tile
 * @sa CHAMELEON_splgtr_Tile
 * @sa CHAMELEON_zplgtr_Tile
 *
 */
int CHAMELEON_zplgtr_Tile( CHAMELEON_Complex64_t bump, cham_uplo_t uplo,
                           CHAM_desc_t *A,
                           unsigned long long int seed )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplgtr_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zplgtr_Tile_Async( bump, uplo, A, seed, sequence, &request );

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
 * @brief Generate a random trapezoidal matrix. Non-blocking equivalent of CHAMELEON_zplgtr_Tile().
 *
 * May return before the computation is finished.
 * Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal to be sure
 *          to have a positive definite matrix if needed.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper trapezoidal or lower trapezoidal:
 *          = ChamUpper: Only the upper trapezoidal part of A is generated;
 *          = ChamLower: Only the lower trapezoidal part of A is generated;
 *          = ChamUpperLower: the full A is generated. Note that in that case
 *          the matrix is symmetric, and it is thus not recommended to call
 *          this function.
 *
 * @param[out] A
 *          On exit, the random trapezoidal matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
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
 * @sa CHAMELEON_zplgtr
 * @sa CHAMELEON_zplgtr_Tile
 * @sa CHAMELEON_cplgtr_Tile_Async
 * @sa CHAMELEON_dplgtr_Tile_Async
 * @sa CHAMELEON_splgtr_Tile_Async
 * @sa CHAMELEON_zplgtr_Tile_Async
 * @sa CHAMELEON_zplgtr_Tile_Async
 *
 */
int CHAMELEON_zplgtr_Tile_Async( CHAMELEON_Complex64_t  bump,
                                 cham_uplo_t            uplo,
                                 CHAM_desc_t           *A,
                                 unsigned long long int seed,
                                 RUNTIME_sequence_t    *sequence,
                                 RUNTIME_request_t     *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplgtr_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zplgtr_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zplgtr_Tile", "NULL request");
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
        chameleon_error("CHAMELEON_zplgtr_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zplgtr_Tile", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return CHAMELEON_SUCCESS;

    chameleon_pzplgsy( bump, uplo, A, A->m, A->i, A->j, seed, sequence, request );

    return CHAMELEON_SUCCESS;
}
