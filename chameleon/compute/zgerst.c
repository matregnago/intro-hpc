/**
 *
 * @file zgerst.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines to restore the original precision of a matrix.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Yuxi Hong
 * @date 2024-02-18
 * @precisions normal z -> d
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Restore the original precision of a given matrix that may have been
 * used in reduced precision during some computations. See
 * CHAMELEON_zgered_Tile() to introduce mixed-precision tiles into the matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of the matrix A:
 *          = ChamUpper: A is upper triangular;
 *          = ChamLower: A is lower triangular;
 *          = ChamUpperLower: A is general.
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix to restore.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgered_Tile
 * @sa CHAMELEON_zgered_Tile_Async
 * @sa CHAMELEON_zgerst_Tile
 * @sa CHAMELEON_cgerst_Tile
 * @sa CHAMELEON_dgerst_Tile
 * @sa CHAMELEON_sgerst_Tile
 *
 */
int CHAMELEON_zgerst_Tile( cham_uplo_t  uplo,
                            CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgerst_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgerst_Tile_Async( uplo, A, sequence, &request );

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
 * @brief Restore the original precision of a given matrix that may have been
 * used in reduced precision during some computations. See
 * CHAMELEON_zgered_Tile() to introduce mixed-precision tiles into the matrix.
 *
 * This is the non-blocking equivalent of CHAMELEON_zgerst_Tile(). It
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
 * @sa CHAMELEON_zgerst_Tile
 * @sa CHAMELEON_zgered_Tile
 * @sa CHAMELEON_zgered_Tile_Async
 *
 */
int CHAMELEON_zgerst_Tile_Async( cham_uplo_t         uplo,
                                  CHAM_desc_t        *A,
                                  RUNTIME_sequence_t *sequence,
                                  RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgerst_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgerst_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgerst_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zgerst_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    chameleon_pzgerst( uplo, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
