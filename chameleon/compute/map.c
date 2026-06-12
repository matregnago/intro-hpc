/**
 *
 * @file map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include "control/common.h"

struct map_args_s {
    cham_unary_operator_t  function;
    void                  *args;
};

static inline int
map_cpu( void *op_args,
         cham_uplo_t uplo, int m, int n, int ndata,
         const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    struct map_args_s *options = (struct map_args_s *)op_args;

    if ( ndata > 1 ) {
        fprintf( stderr, "map_cpu: supports only one piece of data and %d have been given\n", ndata );
    }
    options->function( descA, uplo, m, n, tileA, options->args );

    return 0;
}

static cham_map_operator_t map_op = {
    .name        = "map",
    .cpufunc     = map_cpu,
    .cudafunc    = NULL,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Tile
 *
 *  Apply a given operator on each tile of the given matrix. Operates on
 *  matrices stored by tiles.  All matrices are passed through descriptors.  All
 *  dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] access
 *          - ChamR: A is accessed in read-only mode.
 *          - ChamW: A is accessed in write-only mode.
 *           WARNING: if the descriptor is set for allocation on the fly, the
 *           flush call included in this synchronous API will free all allocated
 *           data, prefer asynchronous call if you want to initialiaze data
 *           before submitting another algorithm.
 *          - ChamRW: A is accessed in read-write mode.
 *
 * @param[in] uplo
 *          - ChamUpper: Only the upper triangular part of the matrix is touched
 *          - ChamLower: Only the lower triangular part of the matrix is touched
 *          - ChamUpperLower: The entire the matrix is touched
 *
 * @param[in,out] A
 *          On exit, the operator has been applied on each tile of the matrix A.
 *
 * @param[in] op_fct
 *          The operator function to apply on each tile of the matrix.
 *
 * @param[in,out] op_args
 *          The arguments structure passed to the operator function when applied
 *          on each tile. May be updated by the operator function.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_map_Tile_Async
 *
 */
int CHAMELEON_map_Tile( cham_access_t         access,
                        cham_uplo_t           uplo,
                        CHAM_desc_t          *A,
                        cham_unary_operator_t op_fct,
                        void                 *op_args )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    {
        cham_map_data_t data = {
            .access = access,
            .desc   = A
        };
        struct map_args_s map_args = {
            .function = op_fct,
            .args     = op_args,
        };
        chameleon_pmap( uplo, 1, &data, &map_op, &map_args, sequence, &request );
    }
    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Tile_Async
 *
 *  Apply a given operator on each tile of the given matrix. Non-blocking equivalent of
 *  CHAMELEON_map_Tile().  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] access
 *          - ChamR: A is accessed in read-only mode.
 *          - ChamW: A is accessed in write-only mode.
 *          INFO: tile of A can be unallocated before the call if the
 *          descriptor is set for allocation on the fly.
 *          - ChamRW: A is accessed in read-write mode.
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
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_map_Tile
 *
 */
int CHAMELEON_map_Tile_Async( cham_access_t         access,
                              cham_uplo_t           uplo,
                              CHAM_desc_t          *A,
                              cham_unary_operator_t op_fct,
                              void                 *op_args,
                              RUNTIME_sequence_t   *sequence,
                              RUNTIME_request_t    *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_map_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_map_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0) {
        return CHAMELEON_SUCCESS;
    }

    {
        cham_map_data_t data = {
            .access = access,
            .desc   = A
        };
        struct map_args_s map_args = {
            .function = op_fct,
            .args     = op_args,
        };
        chameleon_pmap( uplo, 1, &data, &map_op, &map_args, sequence, request );

        /* Need to wait to make sure no one access map_args after this function returned. */
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
    }
    return CHAMELEON_SUCCESS;
}
