/**
 *
 * @file descriptor_pivot.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2025-12-11
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Internal descriptor routines to manipulate LU panel pivot structure.
 *
 */
#include "control/common.h"
#include "control/descriptor.h"
#include "chameleon/runtime.h"

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled descriptor associated to a pivot.
 *
 ******************************************************************************
 *
 * @param[in,out] pivot
 *          The pointer to the pivot descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated pivot descriptor must be generated.
 *
 ******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on success, CHAMELEON_ERR_NOT_INITIALIZED otherwise.
 *
 */
int chameleon_pivot_init( CHAM_desc_pivot_t *pivot, const CHAM_desc_t *desc )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( pivot, 0, sizeof(CHAM_desc_pivot_t) );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    pivot->P    = chameleon_desc_datadist_get_iparam( desc, 0 );
    pivot->Q    = chameleon_desc_datadist_get_iparam( desc, 1 );
    pivot->n    = chameleon_min(desc->mb, desc->nb);
    pivot->nb   = desc->mb;
    pivot->dtyp = desc->dtyp;

    /* Create runtime specific structure like registering data */
    RUNTIME_pivot_create( pivot );

    return rc;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Asynchronously destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] pivot
 *          The pointer to the pivot descriptor to destroy.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated pivot descriptor must be generated.
 *
 */
int chameleon_pivot_destroy_submit( CHAM_desc_pivot_t  *pivot,
                                    RUNTIME_sequence_t *sequence )
{
    RUNTIME_pivot_destroy_submit( sequence, pivot );
    return CHAMELEON_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] pivot
 *          The pointer to the pivot descriptor to destroy.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated pivot descriptor must be generated.
 *
 */
int chameleon_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    RUNTIME_pivot_destroy( pivot );
    return CHAMELEON_SUCCESS;
}
