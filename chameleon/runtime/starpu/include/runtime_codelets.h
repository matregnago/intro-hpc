 /**
 *
 * @file starpu/runtime_codelets.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets main header
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Loris Lucido
 * @author Matthieu Kuhn
 * @date 2025-12-19
 *
 */
#ifndef _runtime_codelets_h_
#define _runtime_codelets_h_

#include "chameleon/config.h"
#include "runtime_codelet_profile.h"

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas.h"

#if defined(CHAMELEON_USE_CUDA)
#include "gpucublas.h"
#endif

#if defined(CHAMELEON_USE_HIP)
#include "gpuhipblas.h"
#endif
#endif /* !defined(CHAMELEON_SIMULATION) */

#if defined(STARPU_CUDA_ASYNC)
#define CODELET_CUDA_FLAGS(_flags_) .cuda_flags = {_flags_},
#else
#define CODELET_CUDA_FLAGS(_flags_)
#endif

#if defined(STARPU_HIP_ASYNC)
#define CODELET_HIP_FLAGS(_flags_) .hip_flags = {_flags_},
#else
#define CODELET_HIP_FLAGS(_flags_)
#endif

#if defined(CHAMELEON_USE_CUDA)
#define CODELET_GPU_FIELDS( gpu_func_name, gpu_flags )                 \
        CODELET_CUDA_FLAGS( gpu_flags )                                \
        .cuda_funcs = {gpu_func_name},
#elif defined(CHAMELEON_USE_HIP)
#define CODELET_GPU_FIELDS( gpu_func_name, gpu_flags )                 \
        CODELET_HIP_FLAGS( gpu_flags )                                 \
        .hip_funcs = {gpu_func_name},
#else
#define CODELET_GPU_FIELDS( gpu_func_name, gpu_flags )
#endif

#define CODELETS_ALL(cl_name, cpu_func_name, gpu_func_name, _original_location_, gpu_flags) \
    struct starpu_perfmodel cl_##cl_name##_model = {                    \
        .type   = STARPU_HISTORY_BASED,                                 \
        .symbol = ""#cl_name                                            \
    };                                                                  \
                                                                        \
    struct starpu_codelet cl_##cl_name = {                              \
        .where     = (_original_location_),                             \
        .cpu_funcs = {(cpu_func_name)},                                 \
        CODELET_GPU_FIELDS( gpu_func_name, gpu_flags )                  \
        .nbuffers  = STARPU_VARIABLE_NBUFFERS,                          \
        .model     = &cl_##cl_name##_model,                             \
        .name      = #cl_name                                           \
    };                                                                  \
                                                                        \
    void cl_##cl_name##_restrict_where(uint32_t where)                  \
    {                                                                   \
        if ( cl_##cl_name.where & where )                               \
            cl_##cl_name.where = (cl_##cl_name.where & where);          \
    }                                                                   \
                                                                        \
    void cl_##cl_name##_restore_where(void)                             \
    {                                                                   \
        cl_##cl_name.where = (_original_location_);                     \
    }

#if defined(CHAMELEON_SIMULATION)
#define CODELETS_CPU(name, cpu_func_name)                    \
    CODELETS_ALL( name, (starpu_cpu_func_t) 1, NULL, STARPU_CPU, 0 )
#else
#define CODELETS_CPU(name, cpu_func_name)                    \
    CODELETS_ALL( name, cpu_func_name, NULL, STARPU_CPU, 0 )
#endif

#if defined(CHAMELEON_USE_HIP)
#define CODELETS_GPU(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_ALL( name, cpu_func_name, gpu_func_name, STARPU_CPU | STARPU_HIP, gpu_flags )
#else
#define CODELETS_GPU(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_ALL( name, cpu_func_name, gpu_func_name, STARPU_CPU | STARPU_CUDA, gpu_flags )
#endif

#define CODELETS_ALL_HEADER(name)                            \
     CHAMELEON_CL_CB_HEADER(name);                           \
     extern struct starpu_codelet cl_##name;                 \
     void cl_##name##_restrict_where(uint32_t where);        \
     void cl_##name##_restore_where(void)

#if defined(CHAMELEON_SIMULATION)

#if defined(CHAMELEON_USE_CUDA) || defined(CHAMELEON_USE_HIP)
#define CODELETS(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_GPU(name, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, gpu_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_CPU(name, (starpu_cpu_func_t) 1)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#endif

#else /* defined(CHAMELEON_SIMULATION) */

#if defined(CHAMELEON_USE_CUDA) //|| defined(CHAMELEON_USE_HIP)
#define CODELETS(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_GPU(name, cpu_func_name, gpu_func_name, gpu_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, cpu_func_name, gpu_func_name, gpu_flags) \
    CODELETS_CPU(name, cpu_func_name)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#endif

#endif

CODELETS_HEADER(map);
CODELETS_HEADER(hgemm);
CODELETS_HEADER(gemmex);

struct cl_hgemm_args_s {
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Real16_t alpha;
    CHAMELEON_Real16_t beta;
};

void
insert_task_convert( const RUNTIME_option_t *options,
                     int m, int n,
                     cham_flttype_t       fromtype,
                     starpu_data_handle_t fromtile,
                     cham_flttype_t       totype,
                     starpu_data_handle_t totile );

#endif /* _runtime_codelets_h_ */
