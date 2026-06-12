/**
 *
 * @file starpu/runtime_workspace.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU workspace header
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _runtime_workspace_h_
#define _runtime_workspace_h_

/*
 * Allocate workspace in host memory: CPU for any worker
 * or allocate workspace in worker's memory: main memory for cpu workers,
 * and embedded memory for CUDA devices.
 */
#define CHAMELEON_HOST_MEM    0
#define CHAMELEON_WORKER_MEM  1

struct chameleon_starpu_ws_s {
    size_t size;
    int    memory_location;
    int    which_workers;
    void  *workspaces[STARPU_NMAXWORKERS];
};

typedef struct chameleon_starpu_ws_s CHAMELEON_starpu_ws_t;

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. CHAMELEON_CUDA|CHAMELEON_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in worker's memory (CHAMELEON_HOST_MEM or CHAMELEON_WORKER_MEM). This function
 * returns 0 upon successful completion.
 */
int   RUNTIME_starpu_ws_alloc   ( CHAMELEON_starpu_ws_t **workspace, size_t size, int which_workers, int memory_location);
int   RUNTIME_starpu_ws_free    ( CHAMELEON_starpu_ws_t  *workspace);
void *RUNTIME_starpu_ws_getlocal( CHAMELEON_starpu_ws_t  *workspace);

#endif /* _runtime_workspace_h_ */
