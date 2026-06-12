/**
 *
 * @file zpotrf_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zpotrf wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_ztile.h"
#if defined(CHAMELEON_USE_CUDA)
#include "gpucublas/gpucublas_z.h"
#endif
#endif

#if !defined(CHAMELEON_SIMULATION)
static inline int
zpotrf_batch_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    cham_uplo_t luplo = (intptr_t)op_args;
    int         info  = 0;
    int         tempnn;

    if ( ndata != 1 ) {
        fprintf( stderr, "zpotrf_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
    }

    tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );

    TCORE_zpotrf(
        luplo, tempnn, tileA, &info );

    if ( info != 0 ) {
        fprintf( stderr, "zpotrf_batch_cpu: Failed to correctly factorize the tile (info = %d)\n", info );
    }

    (void)uplo;
    (void)m;
    return 0;
}
#else
#define zpotrf_batch_cpu NULL
#endif

static cham_map_operator_t zpotrf_batch_map = {
    .name        = "zpotrf",
    .cpufunc     = zpotrf_batch_cpu,
    .cudafunc    = NULL,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zpotrf_batch_Tile - Performs multiple Cholesky factorization in parallel.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A tiles are stored;
 *          = ChamLower: Lower triangle of A tiles are stored.
 *
 * @param[in,out] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *          On exit, each tile is factorized with Cholesky factorization.
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zpotrf_batch_Tile( cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t     data[1];
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrf_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data[0].access = ChamRW;
    data[0].desc   = A;

    chameleon_pmap( ChamUpperLower, 1, data,
                    &zpotrf_batch_map, (void*)((intptr_t)uplo),
                    sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}
