/**
 *
 * @file gpuhipblas.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU kernels main header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @author Loris Lucido
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpuhipblas_h_
#define _gpuhipblas_h_

#include "chameleon/config.h"

#if !defined(CHAMELEON_USE_HIP)
#error "This file should not be included"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>
#include <hip/hip_fp16.h>

#include <hipblas/hipblas.h>
#include <hipsolver/hipsolver.h>

#define HIPBLAS_SADDR(_a_) (&(_a_))
#define HIPBLAS_VALUE(_a_) (_a_)

/**
 * CHAMELEON types and constants
 */
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/constants.h"

/**
 * HIP BLAS headers
 */
BEGIN_C_DECLS

#include "gpuhipblas/gpuhipblas_z.h"
#include "gpuhipblas/gpuhipblas_d.h"
#include "gpuhipblas/gpuhipblas_c.h"
#include "gpuhipblas/gpuhipblas_s.h"

int HIP_hgemm( cham_trans_t transa, cham_trans_t transb,
               int m, int n, int k,
               const hipblasHalf *alpha,
               const hipblasHalf *A, int lda,
               const hipblasHalf *B, int ldb,
               const hipblasHalf *beta,
               hipblasHalf *C, int ldc,
               hipblasHandle_t handle );

END_C_DECLS

/**
 * Coreblas Error
 */
#define gpuhipblas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

extern int chameleon_hipblas_constants[];
#define chameleon_hipblas_const(chameleon_const) chameleon_hipblas_constants[chameleon_const]

END_C_DECLS

#endif /* _gpuhipblas_h_ */
