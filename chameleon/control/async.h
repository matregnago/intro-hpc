/**
 *
 * @file async.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management header
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_async_h_
#define _chameleon_async_h_

#include "chameleon/struct.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
int chameleon_request_fail     (RUNTIME_sequence_t *sequence, RUNTIME_request_t *request, int error);
int chameleon_sequence_create  (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t **sequence);
int chameleon_sequence_destroy (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t *sequence);
int chameleon_sequence_wait    (CHAM_context_t *CHAMELEON, RUNTIME_sequence_t *sequence);
int chameleon_request_create   (CHAM_context_t *CHAMELEON, RUNTIME_request_t **request);
int chameleon_request_destroy  (CHAM_context_t *CHAMELEON, RUNTIME_request_t *request);
int chameleon_request_set      (CHAM_context_t *chamctxt, RUNTIME_request_t *request, int param, int value);

/**
 * @brief Internal function to protect the RUNTIME_data_flush with a test if:
 *    - We are in a recursive task or not
 *    - We are in asynchronous mode and we may not want to flush data for the next
 *      algorithm to run.
 *
 * This wrapper is needed because the flush field is stored in the request and
 * not in the sequence.
 */
static inline void
chameleon_data_flush( const RUNTIME_sequence_t *sequence,
                      const CHAM_desc_t *A, int m, int n, int flush )
{
    if ( !flush ) {
        return;
    }
    RUNTIME_data_flush( sequence, A, m, n );
}

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_async_h_ */
