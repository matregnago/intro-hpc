/**
 *
 * @file descriptor_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users to manipulate IPIV data structures
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "control/descriptor.h"
#include "chameleon/runtime.h"

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] side
 *          Specifies whenever the permutation will be done on the rows or on the columns
 *
 * @param[in] mb
 *          The number of tile in the pivot array.
 *
 * @param[in] m
 *          The size of the pivot array.
 *
 * @param[in] withidx
 *          Boolean to choose if we use the classic perm/invp arrays or the more complex perm/invp_idx arrays that sorts the rows/columns to permute per targeted tile.
 *
 * @param[in] max_m
 *          The maximum value of the ipiv/perm/invp array, i.e. the total number of rows/columns of the matrix to permute. Not referenced if !widthidx.
 *
 * @param[in] max_mt
 *          The maximum number of tiles involved in permutation, i.e. the total number of row/column tiles of the matrix to permute. Not referenced if !widthidx.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] np
 *          The total number of processes.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 * @param[in] get_rankof
 *          The function used to determine which process is responsible for the permutation
 *          of a tile
 ******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on success, CHAMELEON_ERR_NOT_INITIALIZED otherwise.
 *
 */
int chameleon_ipiv_init( CHAM_ipiv_t *ipiv, cham_side_t side, int mb, int m,
                         int withidx, int max_m, int max_mt,
                         int p, int np, void *data,
                         blkrankof_ipiv_fct_t get_rankof )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( ipiv, 0, sizeof(CHAM_ipiv_t) );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ( get_rankof ) {
        ipiv->get_rankof = get_rankof;
    }
    else {
        ipiv->get_rankof = chameleon_getrankof_ipiv_2d_diag;
    }

    ipiv->get_blkdim = chameleon_getblkdim_ipiv;

    ipiv->data    = data;
    ipiv->myrank  = RUNTIME_comm_rank( chamctxt );
    ipiv->i       = 0;
    ipiv->m       = m;
    ipiv->mb      = mb;
    ipiv->mt      = chameleon_ceil( ipiv->m, ipiv->mb );
    ipiv->withidx = withidx;
    ipiv->max_m   = max_m;
    ipiv->max_mt  = max_mt;
    ipiv->P       = p;
    ipiv->NP      = np;

    /* Create runtime specific structure like registering data */
    RUNTIME_ipiv_create( ipiv );

    (void)side;
    return rc;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to destroy.
 *
 */
void chameleon_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    RUNTIME_ipiv_destroy( ipiv );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Create a tiled ipiv descriptor associated to a given matrix.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] side
 *          Specifies whenever the permutation will be done on the rows or on the columns
 *
 * @param[in] desc
 *          The tile descriptor for which an associated ipiv descriptor must be generated.
 *
 * @param[in] m
 *          The size of the pivot array.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS on successful exit
 * @retval CHAMELEON_ERR_NOT_INITIALIZED if failed to initialize the descriptor.
 * @retval CHAMELEON_ERR_OUT_OF_RESOURCES if failed to allocated some ressources.
 *
 */
int CHAMELEON_Ipiv_Create( CHAM_ipiv_t **ipivptr, cham_side_t side, int mb, int m,
                           int p, int np, void *data )
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the ipivriptor */
    ipiv = (CHAM_ipiv_t*)malloc(sizeof(CHAM_ipiv_t));
    if (ipiv == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    chameleon_ipiv_init( ipiv, side, mb, m, 0, -1, -1, p, np, data, NULL );

    *ipivptr = ipiv;
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Destroys an ipiv tile descriptor.
 *
 ******************************************************************************
 *
 * @param[in] ipivptr
 *          The Ipiv tile descriptor to destroy.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Destroy( CHAM_ipiv_t **ipivptr )
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ((ipivptr == NULL) || (*ipivptr == NULL)) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "attempting to destroy a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    ipiv = *ipivptr;
    chameleon_ipiv_destroy( ipiv );
    free(ipiv);
    *ipivptr = NULL;
    return CHAMELEON_SUCCESS;
}

 /**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Flushes the data in the sequence when they won't be reused. This calls
 * cleans up the distributed communication caches, and transfer the data back to
 * the CPU.
 *
 ******************************************************************************
 *
 * @param[in] ipiv
 *          ipiv descriptor.
 *
 * @param[in] sequence
 *          The sequence in which to submit the calls to flush the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Flush( const CHAM_ipiv_t  *ipiv,
                          RUNTIME_sequence_t *sequence )
{
    RUNTIME_ipiv_flushall( sequence, -1, ipiv );
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Gathers an IPIV tile descriptor in a single vector on the given root node.
 *
 ******************************************************************************
 *
 * @param[in] ipivdesc
 *          the ipiv vector descriptor to gather.
 *
 * @param[in] ipiv
 *          The ipiv vector where to store the result. Allocated vector of size
 *          ipivdesc->m on root, not referenced on other nodes.
 *
 * @param[in] root
 *          root node on which to gather the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Gather( const CHAM_ipiv_t *ipivdesc,
                           int               *ipiv,
                           int                root )
{
    CHAM_context_t     *chamctxt = chameleon_context_self();
    RUNTIME_sequence_t *sequence = NULL;
    int                 status;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the tasks to collect the ipiv array on root */
    RUNTIME_ipiv_gather( sequence, ipivdesc, ipiv, root );

    /* Wait for the end */
    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
