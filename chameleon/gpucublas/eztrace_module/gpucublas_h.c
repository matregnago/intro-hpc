/**
 *
 * @file gpucublas_h.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.4.0
 *
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 */
/* -*- c-file-style: "GNU" -*- */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>
#include "gpucublas.h"


static volatile int _gpucublas_h_initialized = 0;

#define CURRENT_MODULE gpucublas_h
DECLARE_CURRENT_MODULE;

#if defined(GPUCUBLAS_HAVE_CUBLASHGEMM)
int (*libCUDA_hgemm)( cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const CHAMELEON_Real16_t *alpha, const CHAMELEON_Real16_t *A, int lda, const CHAMELEON_Real16_t *B, int ldb, const CHAMELEON_Real16_t *beta, CHAMELEON_Real16_t *C, int ldc, cublasHandle_t handle ) = NULL;


int CUDA_hgemm( cham_trans_t transa, cham_trans_t transb,
                int m, int n, int k,
                const CHAMELEON_Real16_t *alpha,
                const CHAMELEON_Real16_t *A, int lda,
                const CHAMELEON_Real16_t *B, int ldb,
                const CHAMELEON_Real16_t *beta,
                CHAMELEON_Real16_t *C, int ldc,
                cublasHandle_t handle ) {
    FUNCTION_ENTRY;
    int ret = libCUDA_hgemm( transa, transb,
                             m, n, k,
                             alpha,
                             A, lda,
                             B, ldb,
                             beta,
                             C, ldc,
                             handle );
    FUNCTION_EXIT;
    return ret;
}
#endif

#if defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
int (*libCUDA_gemmex)( cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const void *alpha, const void *A, int lda, cham_flttype_t Atype, const void *B, int ldb, cham_flttype_t Btype, const void *beta, void *C, int ldc, cham_flttype_t Ctype, cublasHandle_t handle ) = NULL;



int CUDA_gemmex( cham_trans_t transa, cham_trans_t transb,
                 int m, int n, int k,
                 const void *alpha,
                 const void *A, int lda, cham_flttype_t Atype,
                 const void *B, int ldb, cham_flttype_t Btype,
                 const void *beta,
                 void *C, int ldc, cham_flttype_t Ctype,
                 cublasHandle_t handle ) {
    FUNCTION_ENTRY;
    int ret = libCUDA_gemmex( transa, transb,
                              m, n, k,
                              alpha,
                              A, lda, Atype,
                              B, ldb, Btype,
                              beta,
                              C, ldc, Ctype,
                              handle );
    FUNCTION_EXIT;
    return ret;
}
#endif


PPTRACE_START_INTERCEPT_FUNCTIONS(gpucublas_h)
#if defined(CGPUCUBLAS_HAVE_CUBLASHGEMM)
 INTERCEPT3("CUDA_hgemm", libCUDA_hgemm)
#endif
#if defined(CGPUCUBLAS_HAVE_CUBLASGEMMEX)
INTERCEPT3("CUDA_gemmex", libCUDA_gemmex)
#endif
PPTRACE_END_INTERCEPT_FUNCTIONS(gpucublas_h)
