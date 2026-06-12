/**
 *
 * @file struct_context.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context structure
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Xavier Lacoste
 * @date 2025-12-19
 *
 */
#ifndef _struct_context_h_
#define _struct_context_h_

#include "chameleon/types.h"
#include "chameleon/constants.h"
#include "chameleon/runtime_struct.h"

BEGIN_C_DECLS

/**
 *  @brief Chameleon context structure
 */
typedef struct chameleon_context_s {
    RUNTIME_id_t       scheduler;
    int                nworkers;
    int                ncudas;
    int                nthreads_per_worker;

    /* Boolean flags */
    cham_bool_t        warnings_enabled;
    cham_bool_t        autotuning_enabled;
    cham_bool_t        parallel_enabled;
    cham_bool_t        statistics_enabled;
    cham_bool_t        progress_enabled;
    cham_bool_t        generic_enabled;
    cham_bool_t        autominmax_enabled;
    cham_bool_t        optlacpy_enabled;   /**< Enable runtime cpy instead of lacpy kernel            */
    cham_bool_t        runtime_paused;

    cham_householder_t householder;        /**< "domino" (flat) or tree-based (reduction) Householder */
    cham_translation_t translation;        /**< In place or Out of place layout conversion            */

    int                nb;
    int                ib;

    int                first_step;         /**< First step to run in the main algorithm */
    int                last_step;          /**< Last step to run in the main algorithm  */

    int                rhblock;            /**< block size for tree-based (reduction) Householder     */
    int                lookahead;          /**< depth of the look ahead in algorithms                 */
    void              *schedopt;           /**< structure for runtimes                                */
    int                mpi_outer_init;     /**< MPI has been initialized outside our functions        */
    MPI_Comm           comm;               /**< MPI communicator                                      */
} CHAM_context_t;

END_C_DECLS

#endif /* _struct_context_h_ */
