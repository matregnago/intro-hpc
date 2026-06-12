/**
 *
 * @file zgemm_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zgemm wrappers
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

struct zgemm_batch_args_s {
    cham_trans_t          transA;
    cham_trans_t          transB;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};
typedef struct zgemm_batch_args_s zgemm_batch_args_t;

#if !defined(CHAMELEON_SIMULATION)
static inline int
zgemm_batch_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    zgemm_batch_args_t *args = (zgemm_batch_args_t*)op_args;
    const CHAM_desc_t  *descB;
    CHAM_tile_t        *tileB;
    const CHAM_desc_t  *descC;
    CHAM_tile_t        *tileC;
    va_list             ap;
    int                 tempmm, tempnn, tempkk;

    if ( ndata != 3 ) {
        fprintf( stderr, "zgemm_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 3 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descB = va_arg(ap, const CHAM_desc_t *);
    tileB = va_arg(ap, CHAM_tile_t *);
    descC = va_arg(ap, const CHAM_desc_t *);
    tileC = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    tempmm = descC->get_blkdim( descC, m, DIM_m, descC->m );
    tempnn = descC->get_blkdim( descC, n, DIM_n, descC->n );
    if ( args->transA == ChamNoTrans ) {
        tempkk = descA->get_blkdim( descA, n, DIM_n, descA->n );
    }
    else {
        tempkk = descA->get_blkdim( descA, m, DIM_m, descA->m );
    }

    TCORE_zgemm(
        args->transA, args->transB, tempmm, tempnn, tempkk,
        args->alpha, tileA, tileB, args->beta, tileC );

    (void)descB;
    (void)uplo;

    return 0;
}
#else
#define zgemm_batch_cpu NULL
#endif

#if !defined(CHAMELEON_SIMULATION) && defined(CHAMELEON_USE_CUDA)
static inline int
zgemm_batch_cuda( cublasHandle_t handle, void *op_args,
                  cham_uplo_t uplo, int m, int n, int ndata,
                  const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    zgemm_batch_args_t *args = (zgemm_batch_args_t*)op_args;
    const CHAM_desc_t  *descB;
    CHAM_tile_t        *tileB;
    const CHAM_desc_t  *descC;
    CHAM_tile_t        *tileC;
    va_list             ap;
    int                 tempmm, tempnn, tempkk;

    if ( ndata != 3 ) {
        fprintf( stderr, "zgemm_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 3 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descB = va_arg(ap, const CHAM_desc_t *);
    tileB = va_arg(ap, CHAM_tile_t *);
    descC = va_arg(ap, const CHAM_desc_t *);
    tileC = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    tempmm = descC->get_blkdim( descC, m, DIM_m, descC->m );
    tempnn = descC->get_blkdim( descC, n, DIM_n, descC->n );
    if ( args->transA == ChamNoTrans ) {
        tempkk = descA->get_blkdim( descA, n, DIM_n, descA->n );
    }
    else {
        tempkk = descA->get_blkdim( descA, m, DIM_m, descA->m );
    }

    CUDA_zgemm( args->transA, args->transB, tempmm, tempnn, tempkk,
                (cuDoubleComplex*)&(args->alpha),
                tileA->mat, tileA->ld,
                tileB->mat, tileB->ld,
                (cuDoubleComplex*)&(args->beta),
                tileC->mat, tileC->ld,
                handle );

    (void)descB;
    (void)uplo;

    return 0;
}
#else
#define zgemm_batch_cuda NULL
#endif

static cham_map_operator_t zgemm_batch_map = {
    .name        = "zgemm",
    .cpufunc     = zgemm_batch_cpu,
    .cudafunc    = zgemm_batch_cuda,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgemm_batch_Tile - Performs multiple matrix multiplication in parallel.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the tiles from A are transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles from A are not transposed;
 *          = ChamTrans:     tiles from A are transposed;
 *          = ChamConjTrans: tiles from A are conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the tiles from B are transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles from B are not transposed;
 *          = ChamTrans:     tiles from B are transposed;
 *          = ChamConjTrans: tiles from B are conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a collection of mt-by-nt tiles of size C->mb by C->nb
 *          On exit, each tile Cij is overwritten by the matrix:
 *          \f[ alpha * op( A[i,j] )*op( B[i,j] ) * C[i,j] \f]
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zgemm_batch_Tile( cham_trans_t transA, cham_trans_t transB,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                                CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t     data[3];
    zgemm_batch_args_t  params = { transA, transB, alpha, beta };
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data[0].access = ChamR;
    data[0].desc   = A;
    data[1].access = ChamR;
    data[1].desc   = B;
    data[2].access = ( beta == 0. ) ? ChamW : ChamRW;
    data[2].desc   = C;

    chameleon_pmap( ChamUpperLower, 3, data, &zgemm_batch_map, &params, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}
