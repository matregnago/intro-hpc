/**
 *
 * @file starpu/codelet_binit.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Byte init StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2026-03-23
 *
 */
#include "chameleon_starpu_internal.h"

#if !defined(CHAMELEON_SIMULATION)
static void
cl_binit_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t *tileA;
    tileA = cti_interface_get(descr[0]);

    assert( tileA->ld == tileA->m );
    memset( CHAM_tile_get_ptr( tileA ), 0, sizeof(char) * tileA->ld * tileA->n );

    (void)cl_arg;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_binit_cuda_func( void *descr[], void *cl_arg )
{
    cublasStatus_t rc;
    CHAM_tile_t   *tileA;

    tileA = cti_interface_get(descr[0]);
    assert( tileA->ld == tileA->m );

    rc = cudaMemset( CHAM_tile_get_ptr( tileA ), 0, sizeof(char) * tileA->ld * tileA->n );
    assert( rc == CUBLAS_STATUS_SUCCESS );

    (void)cl_arg;
    (void)rc;
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#else /* !defined(CHAMELEON_SIMULATION) */

starpu_cpu_func_t  const cl_binit_cpu_func = (starpu_cpu_func_t)1;
starpu_cuda_func_t const cl_binit_cuda_func = (starpu_cuda_func_t)1;

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
struct starpu_codelet cl_binit = {
#if defined(CHAMELEON_USE_CUDA)
    .where      = STARPU_CPU | STARPU_CUDA,
    .cuda_flags = { STARPU_CUDA_ASYNC },
    .cuda_funcs = { cl_binit_cuda_func },
#else
    .where      = STARPU_CPU,
#endif
    .cpu_funcs  = { cl_binit_cpu_func },
    .nbuffers   = 1,
    .modes      = { STARPU_W },
    .name       = "cl_binit"
};
