/**
 *
 * @file pzbcast.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon broadcast parallel algorithm. These are explicit bcast
 * algorithms when the data is broadcasted to local pieces of data, as opposed
 * to automatic bcast handled by the runtime system where a data A is
 * broadcasted as A.
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A, m, n
#define W(m, n) W, m, n

/**
 * Broadcast the tile A(m,n) to a column buffer with 2d block
 * cyclic distribution
 */
static inline void
__chameleon_pzbcast_tile_bycol_full_2d( const CHAM_desc_t *A,
                                        int                Am,
                                        int                An,
                                        const CHAM_desc_t *W,
                                        int                Wm,
                                        int                Wn,
                                        RUNTIME_option_t  *options )
{
    int tempmm = A->get_blkdim( A, Am, DIM_m, A->m );
    int tempnn = A->get_blkdim( A, An, DIM_n, A->n );
    int P      = chameleon_desc_datadist_get_iparam(W, 0);
    int p;

    for ( p = 0; p < P; p ++ ) {
        INSERT_TASK_zlacpy(
            options,
            ChamUpperLower, tempmm, tempnn,
            A( Am, An ),
            W( Wm + ( ( Am + p ) % P ), Wn ) );
    }
}

/**
 * Broadcast the tile A(m,n) to a column buffer using a ring with 2d block
 * cyclic distribution
 */
static inline void
__chameleon_pzbcast_tile_bycol_ring_2d( const CHAM_desc_t *A,
                                        int                Am,
                                        int                An,
                                        const CHAM_desc_t *W,
                                        int                Wm,
                                        int                Wn,
                                        RUNTIME_option_t  *options )
{
    int tempmm = A->get_blkdim( A, Am, DIM_m, A->m );
    int tempnn = A->get_blkdim( A, An, DIM_n, A->n );
    int P      = chameleon_desc_datadist_get_iparam(W, 0);
    int p;

    INSERT_TASK_zlacpy(
        options,
        ChamUpperLower, tempmm, tempnn,
        A( Am, An ),
        W( Wm + ( Am % P ), Wn ) );

    for ( p = 1; p < P; p ++ ) {
        INSERT_TASK_zlacpy(
            options,
            ChamUpperLower, tempmm, tempnn,
            W( Wm + (( Am + p - 1 ) % P ), Wn ),
            W( Wm + (( Am + p )     % P ), Wn ) );
    }
}

/**
 * Broadcast the tile A(m,n) to a row buffer with 2d block
 * cyclic distribution
 */
static inline void
__chameleon_pzbcast_tile_byrow_full_2d( const CHAM_desc_t *A,
                                        int                Am,
                                        int                An,
                                        const CHAM_desc_t *W,
                                        int                Wm,
                                        int                Wn,
                                        RUNTIME_option_t  *options )
{
    int tempmm = A->get_blkdim( A, Am, DIM_m, A->m );
    int tempnn = A->get_blkdim( A, An, DIM_n, A->n );
    int Q      = chameleon_desc_datadist_get_iparam(W, 1);
    int q;

    for ( q = 0; q < Q; q ++ ) {
        INSERT_TASK_zlacpy(
            options,
            ChamUpperLower, tempmm, tempnn,
            A( Am, An ),
            W( Wm, Wn + (( An + q ) % Q ) ) );
    }
}

/**
 * Broadcast the tile A(m,n) to a row buffer using a ring with 2d block
 * cyclic distribution
 */
static inline void
__chameleon_pzbcast_tile_byrow_ring_2d( const CHAM_desc_t *A,
                                        int                Am,
                                        int                An,
                                        const CHAM_desc_t *W,
                                        int                Wm,
                                        int                Wn,
                                        RUNTIME_option_t  *options )
{
    int tempmm = A->get_blkdim( A, Am, DIM_m, A->m );
    int tempnn = A->get_blkdim( A, An, DIM_n, A->n );
    int Q      = chameleon_desc_datadist_get_iparam(W, 1);
    int q;

    INSERT_TASK_zlacpy(
        options,
        ChamUpperLower, tempmm, tempnn,
        A( Am, An ),
        W( Wm, Wn + ( An % Q ) ) );

    for ( q = 1; q < Q; q ++ ) {
        INSERT_TASK_zlacpy(
            options,
            ChamUpperLower, tempmm, tempnn,
            W( Wm, Wn + (( An + q - 1 ) % Q ) ),
            W( Wm, Wn + (( An + q )     % Q ) ) );
    }
}

typedef void (*chameleon_pzbcast_fct_t)( const CHAM_desc_t *A,
                                         int                Am,
                                         int                An,
                                         const CHAM_desc_t *W,
                                         int                Wm,
                                         int                Wn,
                                         RUNTIME_option_t  *options );

static chameleon_pzbcast_fct_t __chameleon_pzbcast_tile[2][ChamBcastNumber] = {
    /* ChamColumnwise */
    {
        __chameleon_pzbcast_tile_bycol_full_2d,
        __chameleon_pzbcast_tile_bycol_ring_2d,
    },
    /* ChamRowwise */
    {
        __chameleon_pzbcast_tile_byrow_full_2d,
        __chameleon_pzbcast_tile_byrow_ring_2d,
    },
};

/**
 * @brief Bcast a tile A(Am,An) to a workspace W by column or row starting at
 * the position W(Wm,Wn)
 *
 * @param[in] dir
 *          Specify if the tile must be broadcasted on the row or on the column
 *
 * @param[in] A
 *          The descriptor that owns the tile to broadcast
 *
 * @param[in] Am
 *          The row index of the tile to broadcast
 *
 * @param[in] An
 *          The column index of the tile to broadcast
 *
 * @param[in] W
 *          The descriptor of the workspace in which to broadcast the original tile
 *
 * @param[in] Wm
 *          The starting row index of the workspace
 *
 * @param[in] Wn
 *          The starting column index of the workspace
 *
 * @param[in] options
 *          The options argument of the algorithm that needs the broadcast
 *
 */
void
chameleon_pzbcast_tile( cham_store_t       dir,
                        cham_bcast_t       algo,
                        const CHAM_desc_t *A,
                        int                Am,
                        int                An,
                        const CHAM_desc_t *W,
                        int                Wm,
                        int                Wn,
                        RUNTIME_option_t  *options )
{
    assert( A->get_rankof_init == chameleon_getrankof_2d );

    __chameleon_pzbcast_tile[dir-ChamColumnwise][algo]( A, Am, An, W, Wm, Wn, options );
}
