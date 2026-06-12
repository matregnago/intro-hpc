/**
 *
 * @file openmp/runtime_descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP descriptor routines
 *
 * @version 1.4.0
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#include "chameleon_openmp.h"

void *RUNTIME_malloc( size_t size )
{
    return malloc( size );
}

void RUNTIME_free( void  *ptr,
                   size_t size )
{
    (void)size;
    free( ptr );
    return;
}

void RUNTIME_desc_create( CHAM_desc_t *desc )
{
    (void)desc;
    return;
}

void RUNTIME_desc_destroy_submit( CHAM_desc_t *desc, const RUNTIME_sequence_t *sequence )
{
    (void)desc;
    (void)sequence;
    return;
}

void RUNTIME_desc_destroy( CHAM_desc_t *desc )
{
    (void)desc;
    return;
}

int RUNTIME_desc_acquire( const CHAM_desc_t *desc )
{
    (void)desc;
    return CHAMELEON_SUCCESS;
}

int RUNTIME_desc_release( const CHAM_desc_t *desc )
{
    (void)desc;
    return CHAMELEON_SUCCESS;
}

void
RUNTIME_desc_flush( CHAM_desc_t              *desc,
                    const RUNTIME_sequence_t *sequence )
{
    (void)desc;
    (void)sequence;
    return;
}


void
RUNTIME_flush( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

void
RUNTIME_data_flush( const RUNTIME_sequence_t *sequence,
                    const CHAM_desc_t *A, int Am, int An )
{
    (void)sequence;
    (void)A;
    (void)Am;
    (void)An;
    return;
}

void RUNTIME_data_unregister( const RUNTIME_sequence_t *sequence,
                              const CHAM_desc_t *A, int Am, int An )
{
    (void)sequence;
    (void)A;
    (void)Am;
    (void)An;
    return;
}

#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_t *A, int Am, int An, int new_rank )
{
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
}
#endif

void *RUNTIME_data_getaddr( const CHAM_desc_t *desc, int m, int n )
{
    return desc->get_blktile( desc, m, n );
}
