/**
 *
 * @file zhered.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhered wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Ana Hourcau
 * @date 2024-07-17
 * @precisions normal z -> z d
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Computes the Cholesky factorization of a symmetric positive definite
 * or Hermitian positive definite matrix with mixed precision.
 *
 * This is the synchronous version of CHAMELEON_zheredinit_Tile_Async().  It
 * operates on matrices stored by tiles with tiles of potentially different
 * precisions.  All matrices are passed through descriptors.  All dimensions are
 * taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U^H*U or A = L*L^H.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhered
 * @sa CHAMELEON_zhered_Tile_Async
 * @sa CHAMELEON_cpotrfmp_Tile
 * @sa CHAMELEON_dpotrfmp_Tile
 * @sa CHAMELEON_spotrfmp_Tile
 * @sa CHAMELEON_zpotrs_Tile
 *
 */
int CHAMELEON_zhered_Tile( cham_uplo_t uplo, double precision, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zheredinit_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zhered_Tile_Async( uplo, precision, A, sequence, &request );

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
 * @brief Computes the Cholesky factorization of a symmetric positive definite
 * or Hermitian positive definite matrix with mixed precision.
 *
 * This is the non-blocking equivalent of CHAMELEON_zhered_Tile().  It
 * operates on matrices stored by tiles with tiles of potentially different
 * precisions.  All matrices are passed through descriptors.  All dimensions are
 * taken from the descriptors. It may return before the computation is
 * finished. This function allows for pipelining operations at runtime.
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
 * @sa CHAMELEON_zhered
 * @sa CHAMELEON_zhered_Tile
 * @sa CHAMELEON_cpotrfmp_Tile_Async
 * @sa CHAMELEON_dpotrfmp_Tile_Async
 * @sa CHAMELEON_spotrfmp_Tile_Async
 * @sa CHAMELEON_zpotrs_Tile_Async
 *
 */
int CHAMELEON_zhered_Tile_Async( cham_uplo_t uplo, double precision, CHAM_desc_t *A,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhered_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zhered_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zhered_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zhered_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zhered_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /*
     * Quick return
     */
    if ( chameleon_max( A->m, A->n ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    if ( precision < 0. ) {
        char *algostr = chameleon_getenv( "CHAMELEON_GERED_ACC" );
        if ( algostr == NULL ) {
            precision = 1e-12;
        }
        else {
            precision = strtod( algostr, NULL );
        }
    }
    chameleon_pzhered( ChamConjTrans, uplo, precision, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
