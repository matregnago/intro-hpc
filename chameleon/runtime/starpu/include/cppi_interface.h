/**
 *
 * @file starpu/cppi_interface.h
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Header to describe the Chameleon pivot panel interface in StarPU
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @date 2024-08-29
 *
 */
#ifndef _cppi_interface_h_
#define _cppi_interface_h_

#include "chameleon_starpu.h"
#include <starpu_data.h>

extern struct starpu_data_interface_ops cppi_ops;
#define CPPI_INTERFACE_ID cppi_ops.interfaceid

struct cppi_interface_s;
typedef struct cppi_interface_s cppi_interface_t;

/**
 * Chameleon pivot interface
 */
struct cppi_interface_s
{
    CHAM_pivot_t                  pivot;     /**< Copy of the pivot structure                  */
    size_t                        arraysize; /**< Allocated size */
    cham_flttype_t                flttype;   /**< Type of the elements of the matrix           */
    int                           has_diag;  /**< Bool to determine if pivot corresponds to diagonal block of current panel */
    int                           h;         /**< Index of the current column being factorized */
    int                           n;         /**< Number of elements in each row               */
    enum starpu_data_interface_id id;        /**< Identifier of the interface                  */
};

void cppi_interface_init();
void cppi_interface_fini();

CHAM_pivot_t *cppi_handle_get( starpu_data_handle_t handle );
void cppi_register( starpu_data_handle_t *handleptr,
                    cham_flttype_t        flttype,
                    int                   n,
                    int64_t               data_tag,
                    int                   data_rank );

void cl_cppi_redux_cpu_func( void *descr[], void *cl_arg );

#if defined(CHAMELEON_DEBUG_STARPU_CPPI_INTERFACE)
static inline void
cppi_display_dbg( cppi_interface_t *cppi_interface, FILE *f, const char *title )
{
    int i;
    double *diagrow, *pivrow;
    diagrow = cppi_interface->pivot.diagrow;
    pivrow  = cppi_interface->pivot.pivrow;

    fprintf( f, "%sn=%2d, h=%2d, has_diag=%2d, m0=%2d, idx=%2d, interf = %p\n",
             title,
             cppi_interface->n,
             cppi_interface->h,
             cppi_interface->has_diag,
             cppi_interface->pivot.blkm0,
             cppi_interface->pivot.blkidx,
             cppi_interface );

    fprintf(stderr, "Diagonal row: " );
    for( i=0; i<cppi_interface->n; i++) {
        fprintf(stderr, "%e ", diagrow[i] );
    }
    fprintf(stderr, "\n" );
    fprintf(stderr, "Piv      row: " );
    for( i=0; i<cppi_interface->n; i++) {
        fprintf(stderr, "%e ", pivrow[i] );
    }
    fprintf(stderr, "\n" );
}
#else
static inline void
cppi_display_dbg( cppi_interface_t *cppi_interface, FILE *f, const char *title )
{
    (void)cppi_interface;
    (void)f;
    (void)title;
    return;
}
#endif
#endif /* _cppi_interface_h_ */
