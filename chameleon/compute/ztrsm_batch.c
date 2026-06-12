/**
 *
 * @file ztrsm_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch ztrsm wrappers
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

struct ztrsm_batch_args_s {
    cham_side_t           side;
    cham_uplo_t           uplo;
    cham_trans_t          transA;
    cham_diag_t           diag;
    CHAMELEON_Complex64_t alpha;
};
typedef struct ztrsm_batch_args_s ztrsm_batch_args_t;

#if !defined(CHAMELEON_SIMULATION)
static inline int
ztrsm_batch_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    ztrsm_batch_args_t *args = (ztrsm_batch_args_t*)op_args;
    const CHAM_desc_t  *descB;
    CHAM_tile_t        *tileB;
    va_list             ap;
    int                 tempmm, tempnn;

    if ( ndata != 2 ) {
        fprintf( stderr, "ztrsm_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 2 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descB = va_arg(ap, const CHAM_desc_t *);
    tileB = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    tempmm = descB->get_blkdim( descB, m, DIM_m, descB->m );
    tempnn = descB->get_blkdim( descB, n, DIM_n, descB->n );

    TCORE_ztrsm(
        args->side, args->uplo, args->transA, args->diag,
        tempmm, tempnn, args->alpha, tileA, tileB );

    (void)descA;
    (void)descB;
    (void)uplo;

    return 0;
}
#else
#define ztrsm_batch_cpu NULL
#endif

#if !defined(CHAMELEON_SIMULATION) && defined(CHAMELEON_USE_CUDA)
static inline int
ztrsm_batch_cuda( cublasHandle_t handle, void *op_args,
                  cham_uplo_t uplo, int m, int n, int ndata,
                  const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    ztrsm_batch_args_t *args = (ztrsm_batch_args_t*)op_args;
    const CHAM_desc_t  *descB;
    CHAM_tile_t        *tileB;
    va_list             ap;
    int                 tempmm, tempnn;

    if ( ndata != 2 ) {
        fprintf( stderr, "ztrsm_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 2 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descB = va_arg(ap, const CHAM_desc_t *);
    tileB = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    tempmm = descB->get_blkdim( descB, m, DIM_m, descB->m );
    tempnn = descB->get_blkdim( descB, n, DIM_n, descB->n );

    CUDA_ztrsm(
            args->side, args->uplo, args->transA, args->diag,
            tempmm, tempnn, (cuDoubleComplex*)&(args->alpha),
            tileA->mat, tileA->ld,
            tileB->mat, tileB->ld,
            handle );

    (void)descA;
    (void)descB;
    (void)uplo;

    return 0;
}
#else
#define ztrsm_batch_cuda NULL
#endif

static cham_map_operator_t ztrsm_batch_map = {
    .name        = "ztrsm",
    .cpufunc     = ztrsm_batch_cpu,
    .cudafunc    = ztrsm_batch_cuda,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_ztrsm_batch_Tile - Performs multiple triangular solves in parallel.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether tiles of A appears on the left or on the right of tiles of X:
 *          = ChamLeft:  A[i,j] * X[i,j]  = B[i,j]
 *          = ChamRight: X[i,j] * A[i,j]  = B[i,j]
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of tiles of A are stored;
 *          = ChamLower: Lower triangle of tiles of A are stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles of A are transposed;
 *          = ChamTrans:     tiles of A are not transposed;
 *          = ChamConjTrans: tiles of A are conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: tiles of A are non unit;
 *          = ChamUnit:    tiles of A are unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in,out] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *          On exit, each tile B[i,j] is overwritten by X[i,j]
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_ztrsm_batch_Tile( cham_side_t side, cham_uplo_t uplo,
                                cham_trans_t trans, cham_diag_t diag,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t     data[2];
    ztrsm_batch_args_t  params = { side, uplo, trans, diag, alpha };
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_batch_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data[0].access = ChamR;
    data[0].desc   = A;
    data[1].access = ChamRW;
    data[1].desc   = B;

    chameleon_pmap( ChamUpperLower, 2, data, &ztrsm_batch_map, &params, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}
