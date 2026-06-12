/**
 *
 * @file zherk_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zherk wrappers
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

struct zherk_batch_args_s {
    cham_uplo_t  uplo;
    cham_trans_t trans;
    double       alpha;
    double       beta;
};
typedef struct zherk_batch_args_s zherk_batch_args_t;

#if !defined(CHAMELEON_SIMULATION)
static inline int
zherk_batch_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    zherk_batch_args_t *args = (zherk_batch_args_t*)op_args;
    const CHAM_desc_t  *descC;
    CHAM_tile_t        *tileC;
    va_list             ap;
    int                 tempnn, tempkk;

    if ( ndata != 2 ) {
        fprintf( stderr, "zherk_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 2 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descC = va_arg(ap, const CHAM_desc_t *);
    tileC = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    tempnn = descC->get_blkdim( descC, n, DIM_n, descC->n );
    if ( args->trans == ChamNoTrans ) {
        tempkk = descA->get_blkdim( descA, n, DIM_n, descA->n );
    }
    else {
        tempkk = descA->get_blkdim( descA, m, DIM_m, descA->m );
    }

    TCORE_zherk(
        args->uplo, args->trans, tempnn, tempkk,
        args->alpha, tileA, args->beta, tileC );

    (void)uplo;

    return 0;
}
#else
#define zherk_batch_cpu NULL
#endif

#if !defined(CHAMELEON_SIMULATION) && defined(CHAMELEON_USE_CUDA)
static inline int
zherk_batch_cuda( cublasHandle_t handle, void *op_args,
                  cham_uplo_t uplo, int m, int n, int ndata,
                  const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    zherk_batch_args_t *args = (zherk_batch_args_t*)op_args;
    const CHAM_desc_t  *descC;
    CHAM_tile_t        *tileC;
    va_list             ap;
    int                 tempnn, tempkk;

    if ( ndata != 2 ) {
        fprintf( stderr, "zherk_batch_cpu: requires two pieces of data and %d have been given\n", ndata );
        if ( ndata < 2 ) {
            return -1;
        }
    }

    /* Get the second desc */
    va_start(ap, tileA);
    descC = va_arg(ap, const CHAM_desc_t *);
    tileC = va_arg(ap, CHAM_tile_t *);
    va_end(ap);

    tempnn = descC->get_blkdim( descC, n, DIM_n, descC->n );
    if ( args->trans == ChamNoTrans ) {
        tempkk = descA->get_blkdim( descA, n, DIM_n, descA->n );
    }
    else {
        tempkk = descA->get_blkdim( descA, m, DIM_m, descA->m );
    }

    CUDA_zherk( args->uplo, args->trans, tempnn, tempkk,
                &(args->alpha), tileA->mat, tileA->ld,
                &(args->beta),  tileC->mat, tileC->ld,
                handle );

    (void)uplo;

    return 0;
}
#else
#define zherk_batch_cuda NULL
#endif

static cham_map_operator_t zherk_batch_map = {
    .name        = "zherk",
    .cpufunc     = zherk_batch_cpu,
    .cudafunc    = zherk_batch_cuda,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zherk_batch_Tile - Performs multiple rank-k update in parallel.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of C tiles are stored;
 *          = ChamLower: Lower triangle of C tiles are stored.
 *
 * @param[in] trans
 *          Specifies whether the tiles of A are transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles of A are not transposed;
 *          = ChamConjTrans: tiles of A are conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *          On exit, each tile B[i,j] is overwritten by
 *          \f[ B = \alpha [ op( A ) \times conjg( op( A )' )] + \beta B \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zherk_batch_Tile( cham_uplo_t uplo, cham_trans_t trans,
                                double alpha, CHAM_desc_t *A,
                                double beta,  CHAM_desc_t *C )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t     data[2];
    zherk_batch_args_t  params = { uplo, trans, alpha, beta };
    int                 status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zherk_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data[0].access = ChamR;
    data[0].desc   = A;
    data[1].access = ( beta == 0. ) ? ChamW : ChamRW;
    data[1].desc   = C;

    chameleon_pmap( ChamUpperLower, 2, data, &zherk_batch_map, &params, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}
