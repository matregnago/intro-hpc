/**
 *
 * @file starpu/runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU options routines
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"

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
    options->withlacpy = !(chamctxt->optlacpy_enabled);
    options->withcuda  = (chamctxt->ncudas > 0);
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
    int ret = 0;
    if ( worker_size > 0 ) {
        CHAM_tile_t tile = {
            .format = CHAMELEON_TILE_FULLRANK,
            .m      = worker_size,
            .n      = 1,
            .ld     = worker_size,
            .mat    = NULL,
        };
        options->ws_wsize = worker_size;
        starpu_cham_tile_register( (starpu_data_handle_t*)(&(options->ws_worker)),
                                   -1, &tile, ChamByte );
    }
    if ( host_size > 0 ) {
        options->ws_hsize = host_size;
        ret = RUNTIME_starpu_ws_alloc((CHAMELEON_starpu_ws_t**)&(options->ws_host),
                                      host_size, CHAMELEON_CUDA, CHAMELEON_HOST_MEM);
    }
    return ret;
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options )
{
    int ret = 0;
    if ( options->ws_worker != NULL ) {
        starpu_data_unregister_submit((starpu_data_handle_t)(options->ws_worker));
        options->ws_worker = NULL;
    }
    if ( options->ws_host != NULL ) {
        starpu_task_wait_for_all();
        ret = RUNTIME_starpu_ws_free( (CHAMELEON_starpu_ws_t*)(options->ws_host) );
        options->ws_host = NULL;
    }
    return ret;
}
