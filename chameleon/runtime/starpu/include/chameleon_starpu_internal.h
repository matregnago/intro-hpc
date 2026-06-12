/**
 *
 * @file starpu/chameleon_starpu_internal.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU runtime header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @author Samuel Thibault
 * @author Loris Lucido
 * @author Terry Cojean
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_starpu_internal_h_
#define _chameleon_starpu_internal_h_

#include "control/common.h"
#include "chameleon_starpu.h"
#include "chameleon/flops.h"

/* Chameleon interfaces for StarPU */
#include "cham_tile_interface.h"
#include "cppi_interface.h"
#include "cpui_interface.h"

/**
 * @brief StarPU specific Chameleon structure stored in the schedopt field
 */
typedef struct CHAM_context_starpu_s
{
    struct starpu_conf                    starpu_conf; /**< StarPU main configuration structure   */
    struct starpu_parallel_worker_config *pw_config;   /**< StarPU parallel workers configuration */
} CHAM_context_starpu_t;

typedef void (*callback_fct_t)(void *);

/**
 * @brief Convert the Chameleon access enum to the StarPU one
 */
static inline enum starpu_data_access_mode
cham_to_starpu_access( cham_access_t mode ) {
    assert( (enum starpu_data_access_mode)ChamR  == STARPU_R  );
    assert( (enum starpu_data_access_mode)ChamW  == STARPU_W  );
    assert( (enum starpu_data_access_mode)ChamRW == STARPU_RW );
    return (enum starpu_data_access_mode)mode;
}

/**
 * @brief Convert the Chameleon access enum to the StarPU one
 */
static inline cham_access_t
starpu_to_cham_access( enum starpu_data_access_mode mode ) {
    assert( (enum starpu_data_access_mode)ChamR  == STARPU_R  );
    assert( (enum starpu_data_access_mode)ChamW  == STARPU_W  );
    assert( (enum starpu_data_access_mode)ChamRW == STARPU_RW );
    return (cham_access_t)(mode | STARPU_RW);
}

#include "runtime_codelets.h"
#include "runtime_profiling.h"
#include "runtime_codelet_profile.h"
#include "runtime_workspace.h"
#include "runtime_rectasks.h"

#if defined(CHAMELEON_PREC_Z)
extern struct starpu_codelet cl_zinit;
#endif
#if defined(CHAMELEON_PREC_C)
extern struct starpu_codelet cl_cinit;
#endif
#if defined(CHAMELEON_PREC_D)
extern struct starpu_codelet cl_dinit;
#endif
#if defined(CHAMELEON_PREC_S)
extern struct starpu_codelet cl_sinit;
#endif
extern struct starpu_codelet cl_binit;

void *RUNTIME_data_getaddr_withconversion( const RUNTIME_option_t *options,
                                           cham_access_t access, cham_flttype_t flttype,
                                           const CHAM_desc_t *A, int m, int n );

/*
 * MPI Redefinitions
 */
#if defined(CHAMELEON_USE_MPI)

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( options->sequence->comm, (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( options->sequence->comm, (_codelet_), ##__VA_ARGS__ )
#endif

#else

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), ##__VA_ARGS__ )
#endif

#endif

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_shm_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_shm_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), ##__VA_ARGS__ )
#endif

#define INSERT_TASK_COMMON_TASK_PARAMS_NOCB             \
    STARPU_PRIORITY,          options->priority,        \
    STARPU_EXECUTE_ON_WORKER, options->workerid,        \
    STARPU_POSSIBLY_PARALLEL, options->parallel

#define INSERT_TASK_COMMON_TASK_PARAMS( _name_ )        \
    INSERT_TASK_COMMON_TASK_PARAMS_NOCB,                \
    STARPU_CALLBACK,          options->profiling ? cl_##_name_##_callback : NULL

/*
 * Enable codelets names
 */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION > 1))
#define CHAMELEON_CODELETS_HAVE_NAME
#endif

/**
 * MPI tag management
 */
void    chameleon_starpu_tag_init( );
int64_t chameleon_starpu_tag_book( int64_t nbtags );
void    chameleon_starpu_tag_release( int64_t min );

/**
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

void RUNTIME_set_reduction_methods(starpu_data_handle_t handle, cham_flttype_t dtyp);

#include "runtime_mpi.h"
#include "runtime_wontuse.h"

static inline starpu_data_handle_t *
chameleon_starpu_data_gethandle( const CHAM_desc_t *A, int m, int n )
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    return ptrtile;
}

#if defined(CHAMELEON_USE_MPI) && defined(HAVE_STARPU_MPI_CACHED_RECEIVE)
static inline int
chameleon_starpu_data_iscached(const CHAM_desc_t *A, int m, int n)
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    if (!(*ptrtile)) {
        return 0;
    }

    return starpu_mpi_cached_receive(*ptrtile);
}

#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {                 \
        if (chameleon_starpu_data_iscached(A, Am, An)) __chameleon_need_submit = 1; } while(0)

#else

#if defined(CHAMELEON_USE_MPI)
#warning "WAR dependencies need starpu_mpi_cached_receive support from StarPU 1.2.1 or greater"
#endif
#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {} while (0)

#endif

#ifdef CHAMELEON_ENABLE_PRUNING_STATS

#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION \
    int __chameleon_exec = 0; \
    int __chameleon_changed = 0;

#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An) \
    if (chameleon_desc_islocal(A, Am, An)) \
        __chameleon_exec = 1;

#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION \
    RUNTIME_total_tasks++; \
    if (__chameleon_exec) \
        RUNTIME_exec_tasks++; \
    else if (__chameleon_need_submit) \
        RUNTIME_comm_tasks++; \
    else if (__chameleon_changed) \
        RUNTIME_changed_tasks++;

#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank) \
    int __chameleon_myrank; \
    RUNTIME_comm_rank(&__chameleon_myrank); \
    __chameleon_exec = (rank) == __chameleon_myrank; \
    __chameleon_changed = 1; \

#else
#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An)
#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)
#endif

#define RUNTIME_BEGIN_ACCESS_DECLARATION        \
    RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)             \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)            \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)              \
    RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION          \
    RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION;

#define INSERT_TASK_COMMON_PARAMETERS( _name_, _nbuffer_ )      \
    struct starpu_data_descr descrs[_nbuffer_];                 \
    struct starpu_mpi_task_exchange_params params;              \
    struct cl_##_name_##_args_s *clargs  = NULL;                \
    struct starpu_codelet       *cl      = &cl_##_name_;        \
    const char                  *cl_name = #_name_;             \
    int                          nbdata  = 0;

#define INSERT_TASK_COMMON_PARAMETERS_EXTENDED( _name_task_, _name_cl_, _name_arg_, _nbuffer_ ) \
    struct starpu_data_descr descrs[_nbuffer_];                                                 \
    struct starpu_mpi_task_exchange_params params;                                              \
    struct cl_##_name_arg_##_args_s *clargs  = NULL;                                            \
    struct starpu_codelet           *cl      = &cl_##_name_cl_;                                 \
    const char                      *cl_name = #_name_task_;                                    \
    int                              nbdata  = 0;

#define INSERT_TASK_COMMON_PARAMETERS_CLNULL( _name_, _nbuffer_ ) \
    struct starpu_data_descr descrs[_nbuffer_];                   \
    struct starpu_mpi_task_exchange_params params;                \
    struct starpu_codelet           *cl      = NULL;              \
    const char                      *cl_name = #_name_;           \
    int                              nbdata  = 0;

/**
 * This section defines the codelet functions to manage MPI cache and data
 * echanges before and after submitting tasks
 */
#if !defined(CHAMELEON_STARPU_USE_INSERT)

/**
 * @brief Internal function to initialize the StarPU paramas structure.
 *
 * @param[in,out] nbdata
 *          On entry the number of data already registered in descrs. On exist,
 *          the counter is updated if the next handle is registered in the
 *          structure.
 *
 * @param[in,out] descrs
 *          The array of starpu data descriptors (handle + mode). On entry, it
 *          is allcoated to the maximum number of data for the task, and
 *          contains the already registered nbdata handles and their associated
 *          modes. On exit, it is updated with the new handle if needed.
 *
 * @param[in] handle
 *          The data handle
 *
 * @param[in] mode
 *          The access mode
 *
 */
static inline void
starpu_cham_register_descr( int                          *nbdata,
                            struct starpu_data_descr     *descrs,
                            starpu_data_handle_t          handle,
                            enum starpu_data_access_mode  mode )
{
    if ( mode & STARPU_NONE ) {
        return;
    }

    descrs[*nbdata].handle = handle;
    descrs[*nbdata].mode   = mode;

    (*nbdata)++;
    return;
}

#if !defined(CHAMELEON_USE_MPI)

/**
 * @brief Empty data structure to mimic the one provided by StarPU when MPI is enabled
 */
struct starpu_mpi_task_exchange_params {
    int do_execute;
};

static inline void
starpu_cham_exchange_init_params( const RUNTIME_option_t                 *options,
                                  struct starpu_mpi_task_exchange_params *params,
                                  int                                     xrank )
{
    params->do_execute = 1;
    (void)options;
    (void)xrank;
}

static inline void
starpu_cham_exchange_handle_before_execution( const RUNTIME_option_t                 *options,
                                              struct starpu_mpi_task_exchange_params *params,
                                              int                                    *nbdata,
                                              struct starpu_data_descr               *descrs,
                                              starpu_data_handle_t                    handle,
                                              enum starpu_data_access_mode            mode )
{
    starpu_cham_register_descr( nbdata, descrs, handle, mode );

    (void)options;
    (void)params;
    return;
}

static inline void
starpu_cham_exchange_tile_before_execution( const RUNTIME_option_t                 *options,
                                            struct starpu_mpi_task_exchange_params *params,
                                            int                                    *nbdata,
                                            struct starpu_data_descr               *descrs,
                                            const CHAM_desc_t                      *A,
                                            int                                     Am,
                                            int                                     An,
                                            enum starpu_data_access_mode            mode )
{
    starpu_cham_register_descr( nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), mode );

    (void)options;
    (void)params;
    return;
}

#define starpu_cham_task_exchange_data_after_execution( ... ) do {} while(0)

#else

/**
 * @brief Internal function to initialize the StarPU paramas structure.
 *
 * @param[in] options
 *          The runtime options used to set common informations such as
 *          communicator, rank, and priority.
 *
 * @param[in,out] params
 *          On entry, the allocated params structure. On exit the fields of the
 *          structure are initialized.
 *
 * @param[in] xrank
 *          The MPI rank that will execute the task. STARPU_MPI_PER_NODE if all
 *          nodes excute it.
 *
 */
static inline void
starpu_cham_exchange_init_params( const RUNTIME_option_t                 *options,
                                  struct starpu_mpi_task_exchange_params *params,
                                  int                                     xrank )
{
    params->me              = options->sequence->myrank;
    params->xrank           = xrank;
    params->priority        = options->priority;
    params->do_execute      = ( xrank == STARPU_MPI_PER_NODE ) || ( xrank == params->me );
    params->exchange_needed = 0;
}

/**
 * @brief Internal wrapper to starpu_mpi_task_exchange_data_before_execution(),
 * that also perform the cache operation done in the CAHMELEON_ACCESS_X() macros
 * in other runtimes.
 *
 * @param[in] options
 *          The options to parameterize the task
 *
 * @param[in] params
 *          The starpu parameters for the exchange functions. Needs to be
 *          initialized by starpu_cham_init_exchange_param() function.
 *
 * @param[in,out] nbdata
 *          On entry the number of data already registered in descrs. On exist,
 *          the counter is updated if the next handle is registered in the
 *          structure.
 *
 * @param[in,out] descrs
 *          The array of starpu data descriptors (handle + mode). On entry, it
 *          is allcoated to the maximum number of data for the task, and
 *          contains the already registered nbdata handles and their associated
 *          modes. On exit, it is updated with the new handle if needed.
 *
 * @param[in] handle
 *          The data handle
 *
 * @param[in] mode
 *          The access mode
 *
 */
static inline void
starpu_cham_exchange_handle_before_execution( const RUNTIME_option_t                 *options,
                                              struct starpu_mpi_task_exchange_params *params,
                                              int                                    *nbdata,
                                              struct starpu_data_descr               *descrs,
                                              starpu_data_handle_t                    handle,
                                              enum starpu_data_access_mode            mode )
{
    if ( mode & STARPU_NONE ) {
        return;
    }

    starpu_cham_register_descr( nbdata, descrs, handle, mode );

    starpu_mpi_exchange_data_before_execution( options->sequence->comm,
                                               handle, mode, params );

    return;
}

/**
 * @brief Internal wrapper to starpu_mpi_task_exchange_data_before_execution(),
 * that also perform the cache operation done in the CAHMELEON_ACCESS_X() macros
 * in other runtimes.
 *
 * @param[in] options
 *          The options to parameterize the task
 *
 * @param[in] params
 *          The starpu parameters for the exchange functions. Needs to be
 *          initialized by starpu_cham_init_exchange_param() function.
 *
 * @param[in,out] nbdata
 *          On entry the number of data already registered in descrs. On exist,
 *          the counter is updated if the next handle is registered in the
 *          structure.
 *
 * @param[in,out] descrs
 *          The array of starpu data descriptors (handle + mode). On entry, it
 *          is allcoated to the maximum number of data for the task, and
 *          contains the already registered nbdata handles and their associated
 *          modes. On exit, it is updated with the new handle if needed.
 *
 * @param[in] A
 *          The descriptor in which to find the piece of data
 *
 * @param[in] Am
 *          The row index of the piece of data
 *
 * @param[in] An
 *          The column index of the piece of data
 *
 * @param[in] mode
 *          The access mode
 *
 */
static inline void
starpu_cham_exchange_tile_before_execution( const RUNTIME_option_t                 *options,
                                            struct starpu_mpi_task_exchange_params *params,
                                            int                                    *nbdata,
                                            struct starpu_data_descr               *descrs,
                                            const CHAM_desc_t                      *A,
                                            int                                     Am,
                                            int                                     An,
                                            enum starpu_data_access_mode            mode )
{
    unsigned              need_submit = params->do_execute;
    starpu_data_handle_t *ptrtile     = chameleon_starpu_data_gethandle( A, Am, An );

    /*
     * Manage local cache through internal function to avoid the creation of
     * handles if not necessary
     */
    if ( chameleon_desc_islocal( A, Am, An ) ) {
        need_submit = 1;
    }
    else {
        if ( *ptrtile && ( mode & STARPU_W ) &&
             starpu_mpi_cached_receive( *ptrtile ) )
        {
            need_submit = 1;
        }
    }
    if ( options->forcesub && ( mode & STARPU_MPI_REDUX ) )
    {
        need_submit = 1;
    }

    if ( !need_submit ) {
        return;
    }

    /*
     * If we need to submit, let's create the data handle and ask StarPU to perform
     * the necessary communications
     */
    starpu_cham_exchange_handle_before_execution( options, params, nbdata, descrs,
                                                  RTBLKADDR( A, ChamComplexDouble, Am, An ), mode );
    return;
}

/**
 * @brief Internal wrapper to starpu_mpi_task_exchange_data_after_execution().
 *
 * @param[in] options
 *          The options to get the communicator.
 *
 * @param[in] params
 *          The structure that stores a few parameters initialized by
 *          starpu_mpi_task_exchange_data_before_execution().
 *
 * @param[in] nbdata
 *          The size of the descr array.
 *
 * @param[in] descrs
 *          The array of the handle with their associated mode. The array is
 *          initialized in starpu_mpi_task_exchange_data_before_execution().
 *
 */
static inline void
starpu_cham_task_exchange_data_after_execution( const RUNTIME_option_t                *options,
                                                struct starpu_mpi_task_exchange_params params,
                                                int                                    nbdata,
                                                struct starpu_data_descr              *descrs )
{
    starpu_mpi_task_exchange_data_after_execution(
        options->sequence->comm, descrs, nbdata, params );
}

#endif

/**
 * @brief Internal function to initialize the task common parts.
 *
 * @param[in] options
 *          The runtime options used to set common informations such as
 *          communicator, rank, and priority.
 *
 * @param[in,out] task
 *          The task for which to complete the initialization
 *
 * @param[in] nbdata
 *          The size of the descr array.
 *
 * @param[in] descrs
 *          The array of the handle with their associated mode. The array is
 *          initialized in starpu_mpi_task_exchange_data_before_execution().
 *
 * @param[in] callback
 *          The profiling callback function pointe used if profiling is enabled.
 *
 */
static inline void
starpu_cham_task_set_options( const RUNTIME_option_t   *options,
                              struct starpu_task       *task,
                              int                       nbdata,
                              struct starpu_data_descr *descrs,
                              callback_fct_t            callback )
{
    int allocated_buffers = 0;
    int i;

    task->priority = options->priority;

#if defined(CHAMELEON_RUNTIME_SYNC)
    task->synchronous = 1;
#endif

    /* Callback for profiling information */
    if ( options->profiling ) {
        task->callback_func = callback;
    }

    /* Specific worker id */
    if ( options->workerid != -1 ) {
        task->workerid                     = options->workerid;
        task->execute_on_a_specific_worker = 1;
    }

    /* Parallel tasks */
    task->possibly_parallel = options->parallel;

    /* Set the where here */
    // task->where; /* Do restriction here */

    task->nbuffers = nbdata;

    /* Dynamic handles */
    starpu_task_insert_data_make_room( task->cl, task, &allocated_buffers, 0, task->nbuffers );

    for ( i = 0; i < task->nbuffers; i++ ) {
        enum starpu_data_access_mode mode = descrs[i].mode;
        assert( descrs[i].handle );

        STARPU_TASK_SET_HANDLE( task, descrs[i].handle, i );
        STARPU_TASK_SET_MODE( task, mode, i );
    }
}
#endif /* !defined(CHAMELEON_STARPU_USE_INSERT) */

/* Disable STARPU_COMMUTE if CHAMELEON_STARPU_COMMUTE is OFF */
#if !defined(CHAMELEON_STARPU_COMMUTE)
#undef STARPU_COMMUTE
#define STARPU_COMMUTE 0
#endif

#endif /* _chameleon_starpu_internal_h_ */
