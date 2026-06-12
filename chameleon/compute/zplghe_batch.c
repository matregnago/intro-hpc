/**
 *
 * @file zplghe_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zplghe wrappers
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

struct zplghe_batch_args_s {
    double                 bump;
    unsigned long long int seed;
};
typedef struct zplghe_batch_args_s zplghe_batch_args_t;

#if !defined(CHAMELEON_SIMULATION)
static inline int
zplghe_batch_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    zplghe_batch_args_t *args = (zplghe_batch_args_t*)op_args;
    int                 tempnn, m0;

    if ( ndata != 1 ) {
        fprintf( stderr, "zplghe_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
    }

    tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );

    /**
     * Let's fo the math to give the right bigM:
     * jump for the first value is defined as j = m0 + n0 * bigM
     * We need to set (m*n) matrices of size A->mb*A->nb, and we want j, m0, n0 to be defined as:
     *   j = m0 = n0 = (n * A->mt + m) * (A->mb * A->nb)
     * Thus:
     *   bigM = 0;
     */
    m0 = ( n * descA->mt + m ) * (descA->mb * descA->nb );
    TCORE_zplghe( args->bump, tempnn, tempnn, tileA,
                  0, m0, m0, args->seed );

    (void)uplo;

    return 0;
}
#else
#define zplghe_batch_cpu NULL
#endif

static cham_map_operator_t zplghe_batch_map = {
    .name        = "zplghe",
    .cpufunc     = zplghe_batch_cpu,
    .cudafunc    = NULL,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zplghe_batch_Tile - Performs multiple hermitian matrix generation
 *  in parallel.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal of each tile to be sure
 *          they are positive definite matrices.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 * @param[in,out] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *          On exit, each tile is initialized as an hermitian matrix.
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zplghe_batch_Tile( double bump, CHAM_desc_t *A,
                                 unsigned long long int seed )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t     data[1];
    zplghe_batch_args_t params = { bump, seed };
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data[0].access = ChamW;
    data[0].desc   = A;

    chameleon_pmap( ChamUpperLower, 1, data,
                    &zplghe_batch_map, &params,
                    sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}
