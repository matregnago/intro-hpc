/**
 *
 * @file starpu/runtime_wontuse.h
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU wont use implementations to flush pieces of data
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _runtime_wontuse_h_
#define _runtime_wontuse_h_

#include "chameleon_starpu.h"

/**
 * Different implementations of the flush call based on StarPU version
 */
#if defined(HAVE_STARPU_DATA_WONT_USE)

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_wont_use( handle );
}

#elif defined(HAVE_STARPU_IDLE_PREFETCH)

static inline void
chameleon_starpu_data_flush( void *_handle)
{
    starpu_data_handle_t handle = (starpu_data_handle_t)_handle;
    starpu_data_idle_prefetch_on_node( handle, STARPU_MAIN_RAM, 1 );
    starpu_data_release_on_node( handle, -1 );
}

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_acquire_on_node_cb( handle, -1, STARPU_R,
                                    chameleon_starpu_data_flush, handle );
}

#else

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_acquire_cb( handle, STARPU_R,
                            (void (*)(void*))&starpu_data_release, handle );
}

#endif
#endif /* _runtime_wontuse_h_ */
