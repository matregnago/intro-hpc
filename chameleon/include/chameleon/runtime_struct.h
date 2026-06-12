/**
 *
 * @file runtime_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Runtime structures
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Virouleau
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_runtime_struct_h_
#define _chameleon_runtime_struct_h_

#include "chameleon/types.h"

/**
 * Datatype for distributed version
 */
#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#else
#ifndef MPI_Comm
typedef uintptr_t MPI_Comm;
#endif
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif
#endif

BEGIN_C_DECLS

/**
 * @brief Ids of the runtime supported by the RUNTIME API
 */
typedef enum runtime_id_e {
    RUNTIME_SCHED_QUARK,  /**< Quark runtime  */
    RUNTIME_SCHED_PARSEC, /**< PaRSEC runtime */
    RUNTIME_SCHED_STARPU, /**< StarPU runtime */
    RUNTIME_SCHED_OPENMP, /**< OpenMP runtime */
} RUNTIME_id_t;

/**
 * @brief Ids of the worker type
 */
#define RUNTIME_CPU  ((1ULL)<<1)
#define RUNTIME_CUDA ((1ULL)<<2)
#define RUNTIME_HIP  ((1ULL)<<3)

/**
 * @brief Minimal priority value
 */
#define RUNTIME_PRIORITY_MIN 0

/**
 * @brief Maximal priority value
 */
#define RUNTIME_PRIORITY_MAX INT_MAX

/**
 * @brief RUNTIME request structure
 *
 * A request is used to uniquely identifies a set of submitted tasks together,
 * as for example each asynchronous function call.
 *
 */
typedef struct runtime_request_s {
    int   status;   /**< Return status registered by the tasks for the request   */
    int   flush;    /**< Enable/disable automatic data flush within this request */
    int   priority; /**< Field to set a default priority in an algorithm         */
    void *parent;   /**< Specify the parent in recursive submission mode         */
} RUNTIME_request_t;

/**
 *  @brief Runtime request initializer
 */
#define RUNTIME_REQUEST_INITIALIZER { .status = 0, .flush = CHAMELEON_TRUE, .priority = RUNTIME_PRIORITY_MIN, .parent = NULL }

/**
 * @brief RUNTIME sequence structure
 *
 * A sequence is used to uniquely identifies a set of asynchronous function
 * calls sharing common exception handling. If a tasks fails in a request all
 * subsequent tasks in the request, and the sequence will be canceled.
 */
typedef struct runtime_sequence_s {
    int                status;   /**< Return status registered by the tasks for the request     */
    int                myrank;   /**< MPI Comm rank within the associated communicator          */
    RUNTIME_request_t *request;  /**< Pointer to the request that failed if any, NULL otherwise */
    void              *schedopt; /**< Specific runtime data pointer to handle the sequence      */
    MPI_Comm           comm;     /**< MPI communicator                                          */
} RUNTIME_sequence_t;

/**
 * @brief RUNTIME options structure
 *
 * This structure gathers all optionnal fields that can be passed to the runtime
 * system.
 */
typedef struct runtime_option_s {
    RUNTIME_sequence_t *sequence;  /**< Runtime sequence to which attach the submitted tasks     */
    RUNTIME_request_t  *request;   /**< Runtime request to which attach the submitted tasks      */
    int                 profiling; /**< Enable/Disable the profiling of the submitted tasks      */
    int                 parallel;  /**< Enable/Disable the parallel version of submitted tasks   */
    int                 priority;  /**< Define the submitted task priority                       */
    int                 workerid;  /**< Define the prefered worker id to perform the tasks       */
    int                 forcesub;  /**< Force task submission if true                            */
    int                 withlacpy; /**< Flag to force the use of lacpy copies                    */
    int                 withcuda;  /**< Flag to know if cuda is enabled/disabled                 */
    size_t              ws_wsize;  /**< Define the worker workspace size                         */
    size_t              ws_hsize;  /**< Define the host workspace size for hybrid CPU/GPU kernel */
    void               *ws_worker; /**< Pointer to the worker workspace (structure)              */
    void               *ws_host;   /**< Pointer to the host workspace (structure)                */
    void               *schedopt;  /**< Specific runtime data pointer to handle the sequence     */
} RUNTIME_option_t;

END_C_DECLS

#endif /* _chameleon_runtime_struct_h_ */
