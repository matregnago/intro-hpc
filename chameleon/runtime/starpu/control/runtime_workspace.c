/**
 *
 * @file starpu/runtime_workspace.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU workspaces routines
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Loris Lucido
 * @author Philippe Swartvagher
 * @date 2024-09-17
 *
 */
#include "chameleon_starpu_internal.h"

#if defined(CHAMELEON_USE_CUDA)
#define GPU_WORKER_TYPE STARPU_CUDA_WORKER
#define gpuMallocHost( _ptr_, _size_ ) cudaMallocHost( (_ptr_), (_size_) )
#define gpuMalloc       cudaMalloc
#define gpuFreeHost     cudaFreeHost
#define gpuFree         cudaFree
#elif defined(CHAMELEON_USE_HIP)
#define GPU_WORKER_TYPE STARPU_HIP_WORKER
#define gpuMallocHost( _ptr_, _size_ ) hipHostMalloc( (_ptr_), (_size_), 0 )
#define gpuMalloc       hipMalloc
#define gpuFreeHost     hipHostFree
#define gpuFree         hipFree
#endif

static void RUNTIME_allocate_workspace_on_workers(void *arg)
{
    struct chameleon_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;

    int id = starpu_worker_get_id();

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA) || defined(CHAMELEON_USE_HIP)
    type = starpu_worker_get_type(id);
    if ( type == GPU_WORKER_TYPE )
    {
        int memory_location = workspace->memory_location;

        if (memory_location == CHAMELEON_HOST_MEM)
        {
            /* Use pinned memory because the kernel is very likely
             * to transfer these data between the CPU and the GPU.
             * */
            gpuMallocHost( &workspace->workspaces[id], workspace->size );
        }
        else {
            /* Allocate on the device */
            gpuMalloc( &workspace->workspaces[id], workspace->size );
        }
    }
    else
#endif
#endif
    {
        /* This buffer should only be used within the CPU kernel, so
         * there is no point in using pinned memory here. */
        workspace->workspaces[id] = malloc(workspace->size);
    }

    assert( workspace->workspaces[id] );
}


static void RUNTIME_free_workspace_on_workers(void *arg)
{
    struct chameleon_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;
    int id = starpu_worker_get_id();

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA) || defined(CHAMELEON_USE_HIP)
    type = starpu_worker_get_type( id );
    if ( type == GPU_WORKER_TYPE )
    {
        int memory_location = workspace->memory_location;

        if (memory_location == CHAMELEON_HOST_MEM)
        {
            gpuFreeHost( workspace->workspaces[id] );
        }
        else {
            gpuFree( workspace->workspaces[id] );
        }
    }
    else
#endif
#endif
    {
        free( workspace->workspaces[id] );
    }

    workspace->workspaces[id] = NULL;
}

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. CHAMELEON_CUDA|CHAMELEON_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in GPU memory (CHAMELEON_HOST_MEM or CHAMELEON_GPU_MEM). This function
 * returns 0 upon successful completion.:
 */
int RUNTIME_starpu_ws_alloc(CHAMELEON_starpu_ws_t **workspace,
                            size_t size, int which_workers, int memory_location)
{
    if (!workspace)
        return -EINVAL;

    struct chameleon_starpu_ws_s *descr = calloc(1, sizeof(struct chameleon_starpu_ws_s));

    *workspace = descr;

    if (!descr)
        return -ENOMEM;

    descr->size = size;
    descr->memory_location = memory_location;
    descr->which_workers = which_workers;

    starpu_execute_on_each_worker_ex(RUNTIME_allocate_workspace_on_workers, descr, which_workers, "chameleon_ws_alloc");

    return 0;
}

int RUNTIME_starpu_ws_free(CHAMELEON_starpu_ws_t *workspace)
{
    if (!workspace)
        return -EINVAL;

    starpu_execute_on_each_worker_ex(RUNTIME_free_workspace_on_workers, workspace, workspace->which_workers, "chameleon_ws_free");

    free(workspace);

    return 0;
}

void *RUNTIME_starpu_ws_getlocal(CHAMELEON_starpu_ws_t *workspace)
{
    struct chameleon_starpu_ws_s *descr = workspace;
    int id = starpu_worker_get_id();
    return descr->workspaces[id];
}
