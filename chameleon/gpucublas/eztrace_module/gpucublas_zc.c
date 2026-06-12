/**
 *
 * @file gpucublas_zc.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Eztrace module file generated from gpucublas_zc.h
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions mixed zc -> zc ds
 *
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>
#include "gpucublas.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _gpucublas_zc_initialized = 0;

#define CURRENT_MODULE gpucublas_zc
DECLARE_CURRENT_MODULE;

/* pointers to actual gpucublas_zc functions */

#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
int (*libCUDA_clag2z)(int m, int n, const cuFloatComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) = NULL;
int (*libCUDA_zlag2c)(int m, int n, const cuDoubleComplex * A, int lda, cuFloatComplex * B, int ldb, cublasHandle_t handle) = NULL;
#endif



/* Wrapper functions */

#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)

int CUDA_clag2z(int m, int n, const cuFloatComplex * A, int lda, cuDoubleComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_clag2z(m, n, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}

int CUDA_zlag2c(int m, int n, const cuDoubleComplex * A, int lda, cuFloatComplex * B, int ldb, cublasHandle_t handle) {
    FUNCTION_ENTRY;
    int ret = libCUDA_zlag2c(m, n, A, lda, B, ldb, handle);
    FUNCTION_EXIT;
    return ret;
}
#endif


PPTRACE_START_INTERCEPT_FUNCTIONS(gpucublas_zc)
#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
INTERCEPT3("CUDA_clag2z", libCUDA_clag2z)
INTERCEPT3("CUDA_zlag2c", libCUDA_zlag2c)
#endif

PPTRACE_END_INTERCEPT_FUNCTIONS(gpucublas_zc)
