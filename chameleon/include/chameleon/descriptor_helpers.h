/**
 *
 * @file descriptor_helpers.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Set of functions to help the user to declare matrix descriptors (allocation, mapping... )
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @author Lionel Eyraud-Dubois
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 * @addtogroup chameleon_descriptors
 * @{
 *   @brief Set of predefined functions to produce standard matrix mappings
 *
 *   This module provides the set of functions to produce standard mapping of
 *   the matrix tiles among processors, as well as associated function to: get a
 *   tile address within a specific storage, get a tile leading dimension, ...
 *
 */
#ifndef _chameleon_descriptor_helpers_h_
#define _chameleon_descriptor_helpers_h_

#include <chameleon/struct.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name Mapping functions
 * @{
 */
int chameleon_getrankof_2d          ( const CHAM_desc_t *A,    int m, int n );
int chameleon_getrankof_2d_diag     ( const CHAM_desc_t *A,    int m, int n );
int chameleon_getrankof_ipiv_2d_row ( const CHAM_ipiv_t *ipiv, int m, int n );
int chameleon_getrankof_ipiv_2d_col ( const CHAM_ipiv_t *ipiv, int m, int n );
int chameleon_getrankof_ipiv_2d_diag( const CHAM_ipiv_t *ipiv, int m, int n );

typedef struct custom_dist_s{
    int *blocks_dist;         // Matrix of size dist_m times dist_n with values from 1 to number of process MPI
    int dist_m, dist_n;       // The matrix has dist_m rows of dist_n elements
    const char* dist_file;    // Name of the file that contains the distribution
} custom_dist_t;

int chameleon_getrankof_custom_init   ( custom_dist_t **dist, const char *filename );
int chameleon_getrankof_custom_destroy( custom_dist_t **dist );
int chameleon_getrankof_custom        ( const CHAM_desc_t *A, int m, int n );

/**
 * @}
 * @name Panel involvement functions
 * @{
 */

int chameleon_involved_in_panelk_2dbc( const CHAM_desc_t *A, int An );
int chameleon_p_involved_in_panelk_2dbc( const CHAM_desc_t *A, int k, int p );
void chameleon_get_proc_involved_in_panelk_2dbc( const CHAM_desc_t *A,
                                                 int                k,
                                                 int                n,
                                                 void              *ws_reduce );

void chameleon_get_proc_involved_in_rowpanelk_2dbc( const CHAM_desc_t *A,
                                                    int                m,
                                                    int                k,
                                                    void              *ws_reduce );

/**
 * @}
 * @name Block address functions
 * @{
 */
void* chameleon_getaddr_cm  ( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_ccrb( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_null( const CHAM_desc_t *A, int m, int n );
void* chameleon_getaddr_diag( const CHAM_desc_t *A, int m, int n );

/**
 * @}
 * @name Leading dimension functions
 * @{
 */
int chameleon_getblkldd_cm  ( const CHAM_desc_t *A, int m );
int chameleon_getblkldd_ccrb( const CHAM_desc_t *A, int m );
/**
 * @}
 * @name Tile dimensions computation in algorithms
 * @{
 */

/**
 *
 * @ingroup Descriptor
 *
 * @brief Return tile dimension along the m dimension with regular tile sizes
 *
 * @param[in] A
 *          The chameleon descriptor for which to compute the size
 *
 * @param[in] m
 *          The row index of the tile
 *
 * @param[in] lm
 *          The matrix row dimension against which to compute the size
 *
 * @retval The dimension of the tile along the row/first dimension with a limit on lm
 *
 */
static inline int
chameleon_getblkdim_m( const CHAM_desc_t *A, int m, int lm )
{
    return (((m + 1) * A->mb) > lm ) ? lm - m * A->mb : A->mb;
}

/**
 *
 * @ingroup Descriptor
 *
 * @brief Return tile dimension along the n dimension with regular tile sizes
 *
 * @param[in] A
 *          The chameleon descriptor for which to compute the size
 *
 * @param[in] n
 *          The column index of the tile
 *
 * @param[in] ln
 *          The matrix column dimension against which to compute the size
 *
 * @retval The dimension of the tile along the column/second dimension with a limit on ln
 *
 */
static inline int
chameleon_getblkdim_n( const CHAM_desc_t *A, int n, int ln )
{
    return (((n + 1) * A->nb) > ln ) ? ln - n * A->nb : A->nb;
}

/**
 *
 * @ingroup Descriptor
 *
 * @brief Return tile dimension along the dim dimension with regular tile sizes
 *
 * @param[in] A
 *          The chameleon descriptor for which to compute the size
 *
 * @param[in] m
 *          The index of the tile in the given dimension
 *
 * @param[in] dim
 *          The dimension on which to compute the size
 *
 * @param[in] lm
 *          The matrix dimension along the chosen dim.
 *
 * @retval The dimension of the tile along the dim dimension with a limit on lm
 *
 */
static inline int
chameleon_getblkdim( const CHAM_desc_t *A, int m, cham_dim_t dim, int lm )
{
    if ( dim == 0 ) {
        return chameleon_getblkdim_m( A, m, lm );
    }
    else {
        return chameleon_getblkdim_n( A, m, lm );
    }
}

/**
 *
 * @ingroup Descriptor
 *
 * @brief Return tile dimension along the n dimension with regular tile sizes.
 *
 * @param[in] IPIV
 *          The ipiv descriptor.
 *
 * @param[in] m
 *          The index of the tile.
 *
 * @retval The length of the tile.
 *
 */
static inline int
chameleon_getblkdim_ipiv( const CHAM_ipiv_t *IPIV, int m )
{
    return ( ( m + 1 ) * IPIV->mb > IPIV->m ) ? IPIV->m - m * IPIV->mb : IPIV->mb;
}

/**
 * @}
 */
#ifdef __cplusplus
}
#endif

#endif /* _chameleon_descriptor_helpers_h_ */
