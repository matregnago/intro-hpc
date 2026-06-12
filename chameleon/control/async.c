/**
 *
 * @file async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management routines
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 ***
 *
 * @defgroup Sequences
 * @brief Group routines exposed to users to handle asynchronous tasks execution
 *
 */
#include "control/common.h"
#include "chameleon/runtime.h"

/**
 *  Register an exception.
 */
int chameleon_request_fail(RUNTIME_sequence_t *sequence, RUNTIME_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return status;
}

/**
 *  Create a sequence
 */
int chameleon_sequence_create(CHAM_context_t *chamctxt, RUNTIME_sequence_t **sequence)
{
    if ((*sequence = malloc(sizeof(RUNTIME_sequence_t))) == NULL) {
        chameleon_error("CHAMELEON_Sequence_Create", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    (*sequence)->comm   = chamctxt->comm;
    (*sequence)->myrank = RUNTIME_comm_rank( chamctxt );

    RUNTIME_sequence_create( chamctxt, *sequence );

    (*sequence)->status = CHAMELEON_SUCCESS;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int chameleon_sequence_destroy(CHAM_context_t *chamctxt, RUNTIME_sequence_t *sequence)
{
    RUNTIME_sequence_destroy( chamctxt, sequence );
    free(sequence);
    return CHAMELEON_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int chameleon_sequence_wait(CHAM_context_t *chamctxt, RUNTIME_sequence_t *sequence)
{
    if ( chamctxt->runtime_paused ) {
        chameleon_warning("chameleon_sequence_wait", "CHAMELEON was paused, so it is resumed and it will be paused after the wait");
        RUNTIME_resume( chamctxt );
    }
    RUNTIME_sequence_wait( chamctxt, sequence );
    if ( chamctxt->runtime_paused ) {
        RUNTIME_pause( chamctxt );
    }
    assert( sequence->status == CHAMELEON_SUCCESS );
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Sequences
 *
 *  CHAMELEON_Sequence_Create - Create a squence.
 *
 ******************************************************************************
 *
 * @param[out] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Sequence_Create(RUNTIME_sequence_t **sequence)
{
    CHAM_context_t *chamctxt;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    status = chameleon_sequence_create(chamctxt, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  CHAMELEON_Sequence_Destroy - Destroy a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Sequence_Destroy(RUNTIME_sequence_t *sequence)
{
    CHAM_context_t *chamctxt;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Destroy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Destroy", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    status = chameleon_sequence_destroy(chamctxt, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  CHAMELEON_Sequence_Wait - Wait for the completion of a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Sequence_Wait(RUNTIME_sequence_t *sequence)
{
    CHAM_context_t *chamctxt;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Wait", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Wait", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    status = chameleon_sequence_wait(chamctxt, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  CHAMELEON_Sequence_Flush - Terminate a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 * @param[in] request
 *          The flush request.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Sequence_Flush(RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Flush", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_Sequence_Flush", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    RUNTIME_sequence_flush( chamctxt->schedopt, sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);

    return CHAMELEON_SUCCESS;
}
