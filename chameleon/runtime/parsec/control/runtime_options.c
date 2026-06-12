/**
 *
 * @file parsec/runtime_options.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC options routines
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include "chameleon_parsec.h"

void RUNTIME_options_init( RUNTIME_option_t     *options,
                           const CHAM_context_t *chamctxt,
                           RUNTIME_sequence_t   *sequence,
                           RUNTIME_request_t    *request )
{
    options->sequence  = sequence;
    options->request   = request;
    options->profiling = CHAMELEON_STATISTICS == CHAMELEON_TRUE;
    options->parallel  = CHAMELEON_PARALLEL == CHAMELEON_TRUE;
    options->priority  = request->priority;
    options->workerid  = -1;
    options->forcesub  = 0;
    options->withlacpy = 0;
    options->ws_wsize  = 0;
    options->ws_hsize  = 0;
    options->ws_worker = NULL;
    options->ws_host   = NULL;
    return;
}

void RUNTIME_options_finalize( RUNTIME_option_t     *options,
                               const CHAM_context_t *chamctxt )
{
    (void)options;
    (void)chamctxt;
    return;
}

int RUNTIME_options_ws_alloc( RUNTIME_option_t *options, size_t worker_size, size_t host_size )
{
    options->ws_wsize = worker_size;
    options->ws_hsize = host_size;
    return CHAMELEON_SUCCESS;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    options->ws_wsize = 0;
    options->ws_hsize = 0;
    return CHAMELEON_SUCCESS;
}
