/**
 *
 * @file control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon control routines
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Virouleau
 * @author Samuel Thibault
 * @author Philippe Swartvagher
 * @author Loris Lucido
 * @date 2024-03-16
 *
 ***
 *
 * @defgroup Control
 * @brief Group routines exposed to users to control CHAMELEON state
 *
 */
#include "control/common.h"
#include "chameleon/runtime.h"

/**
 *
 * @ingroup Control
 *
 * @brief Initialize CHAMELEON with number of cpus and gpus (using
 * MPI_COMM_WORLD).
 *
 ******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use.
 *
 * @param[in] gpus
 *          Number of cuda devices to use.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int __chameleon_init(int cores, int gpus)
{
   return __chameleon_initpar(cores, gpus, -1);
}

/**
 *
 * @ingroup Control
 *
 * @brief Initialize CHAMELEON with number of cpus and gpus and threads per
 * worker (using MPI_COMM_WORLD).
 *
 ******************************************************************************
 *
 * @param[in] ncpus
 *          Number of cores to use.
 *
 * @param[in] ngpus
 *          Number of cuda devices to use.
 *
 * @param[in] nthreads_per_worker
 *          Number of threads per worker (cpu, cuda device).
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int __chameleon_initpar(int ncpus, int ngpus, int nthreads_per_worker)
{
    return __chameleon_initparcomm( ncpus, ngpus, nthreads_per_worker, MPI_COMM_WORLD );
}

/**
 *
 * @ingroup Control
 *
 * @brief Initialize CHAMELEON with number of cpus and gpus and threads per
 * worker and using a given MPI communicator.
 *
 ******************************************************************************
 *
 * @param[in] ncpus
 *          Number of cores to use.
 *
 * @param[in] ngpus
 *          Number of cuda devices to use.
 *
 * @param[in] nthreads_per_worker
 *          Number of threads per worker (cpu, cuda device).
 *
 * @param[in] comm
 *          The MPI communicator.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int __chameleon_initparcomm(int ncpus, int ngpus, int nthreads_per_worker, MPI_Comm comm)
{
    CHAM_context_t *chamctxt;

    /* Create context and insert in the context map */
    chamctxt = chameleon_context_create();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Init", "chameleon_context_create() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

#if defined(CHAMELEON_USE_MPI)
#  if defined(CHAMELEON_SIMULATION)
    /* Assuming that we don't initialize MPI ourself (which SMPI doesn't support anyway) */
    chamctxt->mpi_outer_init = 1;
#  else
    {
        int flag = 0, provided = 0;
        MPI_Initialized( &flag );
        chamctxt->mpi_outer_init = flag;
        if ( !flag ) {
            /* MPI_THREAD_SERIALIZED should be enough.
             * In testings, only StarPU's internal thread performs
             * communications, and *then* Chameleon performs communications in
             * the check step. */
            const int required = MPI_THREAD_MULTIPLE;
            if ( MPI_Init_thread( NULL, NULL, required, &provided ) != MPI_SUCCESS) {
                chameleon_fatal_error("CHAMELEON_Init", "MPI_Init_thread() failed");
            }
            if ( provided < required ) {
                chameleon_fatal_error("CHAMELEON_Init",
                                      "MPI_Init_thread() was not able to provide the requested thread support (MPI_THREAD_MULTIPLE),\n"
                                      "this may be an issue if the level provided is not enough for the underlying runtime system." );
            }
        }
    }
#  endif
#endif

#if !defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_USE_HIP)
    if ( ngpus != 0 ) {
        chameleon_warning("CHAMELEON_Init", "CHAMELEON_USE_CUDA or CHAMELEON_USE_HIP are not defined, ngpus is forced to 0");
        ngpus = 0;
    }
#endif

    chamctxt->ncudas = ngpus;
    chamctxt->comm = comm;
    return RUNTIME_init( chamctxt, ncpus, ngpus, nthreads_per_worker );
}

/**
 *
 * @ingroup Control
 *
 * @brief Finalize CHAMELEON.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int __chameleon_finalize(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Finalize", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Make sure all data are flushed */
    RUNTIME_flush( chamctxt );

    /* Wait for anything running */
#  if !defined(CHAMELEON_SIMULATION)
    RUNTIME_barrier(chamctxt);
#  endif

    /* Stop the runtime system */
    RUNTIME_finalize( chamctxt );

#if defined(CHAMELEON_USE_MPI)
    /* Finalize MPI if initialized by Chameleon */
    if ( !chamctxt->mpi_outer_init ) {
        MPI_Finalize();
    }
#endif

    chameleon_context_destroy();
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 * @brief Check if the CHAMELEON library is initialized or not.
 *
 * @retval True if initialized
 * @retval False if not initialized
 *
 */
int CHAMELEON_Initialized(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return 0;
    }
    else {
        return 1;
    }
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Pause - Suspend CHAMELEON runtime to poll for new tasks.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Pause(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Pause", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( chamctxt->runtime_paused ) {
        chameleon_warning("CHAMELEON_Pause", "CHAMELEON already paused");
        return CHAMELEON_SUCCESS;
    }
    chamctxt->runtime_paused = CHAMELEON_TRUE;
    RUNTIME_pause(chamctxt);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Resume - Symmetrical call to CHAMELEON_Pause,
 *  used to resume the workers polling for new tasks.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Resume(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Resume", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( !chamctxt->runtime_paused ) {
        chameleon_warning("CHAMELEON_Resume", "CHAMELEON was already resumed");
        return CHAMELEON_SUCCESS;
    }
    chamctxt->runtime_paused = CHAMELEON_FALSE;
    RUNTIME_resume(chamctxt);
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Distributed_start - Prepare the distributed processes for computation
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Distributed_start(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Finalize", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_barrier( chamctxt );
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Distributed_stop - Clean the distributed processes after computation
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Distributed_stop(void)
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Finalize", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    RUNTIME_barrier( chamctxt );
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Comm_size - Return the size of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The size of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int CHAMELEON_Comm_size()
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Comm_size", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_comm_size( chamctxt );
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_Comm_rank - Return the rank of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The rank of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int CHAMELEON_Comm_rank()
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Comm_rank", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_comm_rank( chamctxt );
}

/**
 *
 * @ingroup Control
 *
 *  CHAMELEON_GetThreadNbr - Return the number of CPU workers initialized by the
 *  runtime
 *
 ******************************************************************************
 *
 * @retval The number of CPU workers started
 *
 */
int CHAMELEON_GetThreadNbr( )
{
    CHAM_context_t *chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_GetThreadNbr", "CHAMELEON not initialized");
        return -1;
    }

    return RUNTIME_thread_size( chamctxt );
}
