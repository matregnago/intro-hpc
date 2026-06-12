/**
 *
 * @file starpu/codelet_zinit.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zinit StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2026-03-23
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zinit_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlaset( ChamUpperLower, tileA->m, tileA->n, 0., 0., tileA );

    (void)cl_arg;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zinit_cuda_func( void *descr[], void *cl_arg )
{
    cublasHandle_t  handle = starpu_cublas_get_local_handle();
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zzero = make_cuDoubleComplex(0.0, 0.0);
#else
    double          zzero = 0.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
    CHAM_tile_t    *tileA;

    tileA = cti_interface_get(descr[0]);

    TCUDA_zlaset( ChamUpperLower, tileA->m, tileA->n,
                  &zzero, &zzero,
                  tileA, handle );
    (void)cl_arg;
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#else /* !defined(CHAMELEON_SIMULATION) */

starpu_cpu_func_t const cl_zinit_cpu_func  = (starpu_cpu_func_t)1;
starpu_cuda_func_t const cl_zinit_cuda_func = (starpu_cuda_func_t)1;

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
struct starpu_codelet cl_zinit = {
#if defined(CHAMELEON_USE_CUDA)
    .where      = STARPU_CPU | STARPU_CUDA,
    .cuda_flags = { STARPU_CUDA_ASYNC },
    .cuda_funcs = { cl_zinit_cuda_func },
#else
    .where      = STARPU_CPU,
#endif
    .cpu_funcs  = { cl_zinit_cpu_func },
    .nbuffers   = 1,
    .modes      = { STARPU_W },
    .name       = "cl_zinit"
};
