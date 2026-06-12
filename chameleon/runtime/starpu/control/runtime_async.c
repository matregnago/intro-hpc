/**
 *
 * @file starpu/runtime_async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU asynchronous routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"

/**
 *  Create a sequence
 */
int RUNTIME_sequence_create( CHAM_context_t     *chamctxt,
                             RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int RUNTIME_sequence_destroy( CHAM_context_t     *chamctxt,
                              RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int RUNTIME_sequence_wait( CHAM_context_t     *chamctxt,
                           RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;

    if (chamctxt->progress_enabled) {
        RUNTIME_progress(chamctxt);
    }

#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_WAIT_FOR_ALL)
    starpu_mpi_wait_for_all(sequence->comm);
#  else
    starpu_task_wait_for_all();
    starpu_mpi_barrier(sequence->comm);
#  endif
#else
    starpu_task_wait_for_all();
#endif
    return CHAMELEON_SUCCESS;
}

/**
 *  Terminate a sequence
 */
void RUNTIME_sequence_flush( CHAM_context_t     *chamctxt,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request,
                             int status )
{
    (void)chamctxt;
    assert( status == CHAMELEON_SUCCESS );
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}
