/**
 *
 * @file starpu/cpui_interface.h
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Header to describe the Chameleon permutation interface with StarPU used in laswp
 *
 * @version 1.4.0
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#ifndef _cpui_interface_h_
#define _cpui_interface_h_

#include "chameleon_starpu.h"
#include <starpu_data.h>

extern struct starpu_data_interface_ops cpui_ops;
#define CPUI_INTERFACE_ID STARPU_UNKNOWN_INTERFACE_ID

struct cpui_interface_s;
typedef struct cpui_interface_s cpui_interface_t;

/**
 * Chameleon pivot interface
 */
struct cpui_interface_s
{
    CHAM_laswpws_t                ws;      /**< Copy of the workspace structure                                  */
    cham_flttype_t                flttype; /**< Type of the elements of the matrix                               */
    cham_side_t                   side;    /**< Specifies whether the permutation is done on the rows or columns */
    int                           m;       /**< Number of rows per tile                                          */
    int                           n;       /**< Number of columns per tile                                       */
    enum starpu_data_interface_id id;      /**< Identifier of the interface                                      */
};

void cpui_interface_init();
void cpui_interface_fini();

CHAM_laswpws_t *cpui_handle_get( starpu_data_handle_t handle );
void cpui_register( starpu_data_handle_t *handleptr,
                    cham_side_t           side,
                    cham_flttype_t        flttype,
                    int                   m,
                    int                   n,
                    int64_t               data_tag,
                    int                   data_rank );

void cl_cpui_redux_cpu_func( void *descr[], void *cl_arg );

#if defined(CHAMELEON_DEBUG_STARPU_CPUI_INTERFACE)
static inline void
cpui_display_dbg( cpui_interface_t *cpui_interface, FILE *f, const char *title )
{
    int i;
    int *index = cpui_interface->ws.index;

    fprintf( f, "%sside=%2d, n=%2d, nindex=%2d, interf = %p\n",
             title,
             cpui_interface->n,
             cpui_interface->h,
             cpui_interface->ws.nindex,
             cpui_interface );

    fprintf(stderr, "index: " );
    for( i=0; i<cpui_interface->ws.nindex; i++) {
        fprintf(stderr, "%d ", index[i] );
    }
    fprintf(stderr, "\n" );
}
#else
static inline void
cpui_display_dbg( cpui_interface_t *cpui_interface, FILE *f, const char *title )
{
    (void)cpui_interface;
    (void)f;
    (void)title;
    return;
}
#endif
#endif /* _cpui_interface_h_ */
