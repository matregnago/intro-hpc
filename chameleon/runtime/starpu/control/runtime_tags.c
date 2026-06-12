/**
 *
 * @file starpu/runtime_tags.c
 *
 * @copyright 2017-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.4.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 * Functions to manage the MPI data tags with StarPU (originated from PaStiX).
 *
 * @addtogroup chameleon_starpu
 * @{
 *
 **/
#include "chameleon_starpu_internal.h"

#if defined(CHAMELEON_USE_MPI)
#if defined(HAVE_STARPU_MPI_TAGS_ALLOCATE)
void
chameleon_starpu_tag_init( void )
{
    return;
}

int64_t
chameleon_starpu_tag_book( int64_t nbtags )
{
    return starpu_mpi_tags_allocate( nbtags );
}

void
chameleon_starpu_tag_release( int64_t min )
{
    starpu_mpi_tags_free( min );
}

#else
/**
 * @brief Structure Chameleon StarPU tag
 *
 * List structure to manage the set of available tags.
 */
struct cst_range_;
typedef struct cst_range_ cst_range_t;

struct cst_range_ {
    int64_t      min;  /**< Minimal value in the range     */
    int64_t      max;  /**< Maximal value in the range     */
    cst_range_t *next; /**< Pointer to the following range */
};

/**
 * @brief Pointer to the first set or registered tags
 */
static cst_range_t *cst_first     = NULL;

/**
 * @brief StarPU tag upper bound
 */
static int64_t starpu_tag_ub = 0;

/**
 *******************************************************************************
 *
 * @brief Initialize the StarPU tags manager.
 *
 ******************************************************************************/
void
chameleon_starpu_tag_init( void )
{
    if (!starpu_tag_ub) {
        int          ok       = 0;
        void        *tag_ub_p = NULL;

        CHAM_context_t *chamctxt = chameleon_context_self();
        starpu_mpi_comm_get_attr( chamctxt->comm, STARPU_MPI_TAG_UB, &tag_ub_p, &ok );
        starpu_tag_ub = (uint64_t)((intptr_t)tag_ub_p);

        if ( !ok ) {
            chameleon_error("chameleon_starpu_tag_init", "MPI_TAG_UB not known by StarPU\n");
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Book a range of StarPU unique tags of size nbtags.
 *
 * This function returns the minimal tag value available to allow the
 * registration of nbtags data in a continuous range.
 *
 * Note that this function must be called exactly the same way on all nodes to
 * make sure the tags are identical from one node to another.
 *
 *******************************************************************************
 *
 * @param[in] nbtags
 *          The number of tags required to register the sparse matrix or right
 *          hand side.
 *
 *******************************************************************************
 *
 * @return V, the minimal tag value to use. The range [V:V+nbtags-1] is booked.
 *
 ********************************************************************************/
int64_t
chameleon_starpu_tag_book( int64_t nbtags )
{
    cst_range_t *new;
    cst_range_t *prev    = NULL;
    cst_range_t *current = cst_first;
    int64_t      min = 0;
    int64_t      max = ( current == NULL ) ? starpu_tag_ub : current->min;

    if ( nbtags == 0 ) {
        return -1;
    }
    assert( starpu_tag_ub != 0 ); /* StarPU tag must be initialized */

    while ( ((max - min) < nbtags) && (current != NULL) ) {
        min     = current->max;
        prev    = current;
        current = current->next;
        max     = ( current == NULL ) ? starpu_tag_ub : current->min;
    }

    if ( (max - min) < nbtags ) {
        /* chameleon_error( "chameleon_starpu_tag_book: No space left in tags (looking for %ld tags)\n", */
        /*                  nbtags ); */
        return -1;
    }

    new = malloc( sizeof( cst_range_t ) );
    new->min  = min;
    new->max  = min + nbtags;
    new->next = current;
    if ( prev == NULL ) {
        cst_first = new;
    }
    else {
        assert( prev->next == current );
        prev->next = new;
    }

#if defined(CHAMELEON_DEBUG_STARPU)
    fprintf( stderr, "[%02d] chameleon_starpu_tag: Book %ld - %ld\n",
             CHAMELEON_Comm_rank(), min, min + nbtags );
#endif

    assert( cst_first != NULL );
    return new->min;
}

/**
 *******************************************************************************
 *
 * @brief Release the set of tags starting by min.
 *
 * This function releases the range of tags that starts by the min value.
 *
 *******************************************************************************
 *
 * @param[in] min
 *          The initial value in the range
 *
 ******************************************************************************/
void
chameleon_starpu_tag_release( int64_t min )
{
    cst_range_t *prev    = NULL;
    cst_range_t *current = cst_first;

    assert( cst_first != NULL ); /* At least one range must be registered */
    if ( current == NULL ) {
#if defined(CHAMELEON_DEBUG_STARPU)
        fprintf( stderr, "chameleon_starpu_tag: FAILED to release [%ld,...]\n",
                 min );
#endif

        return;
    }

    while ( (current != NULL) && (current->min < min) ) {
        prev    = current;
        current = current->next;
    }

    if ( current == NULL ) {
#if defined(CHAMELEON_DEBUG_STARPU)
        fprintf( stderr, "chameleon_starpu_tag: FAILED to release [%ld,...] no set registered with this min value\n",
                 min );
#endif

        return;
    }

    assert( current != NULL );
    assert( current->min == min );

    if ( prev ) {
        prev->next = current->next;
    }
    else {
        assert( current == cst_first );
        cst_first = current->next;
    }

#if defined(CHAMELEON_DEBUG_STARPU)
    fprintf( stderr, "chameleon_starpu_tag: Release %ld - %ld\n",
             current->min, current->max );
#endif

    free( current );

    return;
}

#endif /* defined(HAVE_STARPU_MPI_TAGS_ALLOCATE) */
#else /* defined(CHAMELEON_USE_MPI) */

/**
 *******************************************************************************
 *
 * @brief Initialize the StarPU tags manager.
 *
 ******************************************************************************/
void
chameleon_starpu_tag_init( ) {
    return;
}

/**
 *******************************************************************************
 *
 * @brief Book a range of StarPU unique tags of size nbtags.
 *
 * This function returns the minimal tag value available to allow the
 * registration of nbtags data in a continuous range.
 *
 * Note that this function must be called exactly the same way on all nodes to
 * make sure the tags are identical from one node to another.
 *
 *******************************************************************************
 *
 * @param[in] nbtags
 *          The number of tags required to register the sparse matrix or right
 *          hand side.
 *
 *******************************************************************************
 *
 * @return V, the minimal tag value to use. The range [V:V+nbtags-1] is booked.
 *
 ********************************************************************************/
int64_t
chameleon_starpu_tag_book( __attribute__((unused)) int64_t nbtags ) {
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Release the set of tags starting by min.
 *
 * This function releases the range of tags that starts by the min value.
 *
 *******************************************************************************
 *
 * @param[in] min
 *          The initial value in the range
 *
 ******************************************************************************/
void
chameleon_starpu_tag_release( __attribute__((unused)) int64_t min ) {
    return;
}

#endif

/**
 * @}
 */
