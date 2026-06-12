/**
 *
 * @file cuda_dlag2h.cu
 *
 * @copyright 2023-2023 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_dlag2h GPU kernel
 *
 * @version 1.4.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal d -> d s
 *
 * This file is an adaptation of the MAGMA zlag2c.cu, clag2z.cu, hlag2s.cu, and slag2h.cu
 *
 */
#include "gpucublas.h"
#include <coreblas/lapacke.h>

#if CUDA_VERSION < 7500
#error "This file should not be included as the half precision floats are not supported."
#endif

#define BLK_X 64
#define BLK_Y 32

__device__  int cuda_dlag2h_flag = 0;

/*
 * Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
 * Each block has BLK_X threads.
 * Each thread loops across one row, updating BLK_Y entries.
 */
__global__
void cuda_dlag2h_kernel(
    int m, int n,
    const double *A, int lda,
    CHAMELEON_Real16_t *HA, int ldha,
    double rmax )
{
    int ind = blockIdx.x * BLK_X + threadIdx.x;
    int iby = blockIdx.y * BLK_Y;

    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    double tmp;

    /* do only rows inside matrix */
    if ( ind > m ) {
        return;
    }

    A  += ind + iby*lda;
    HA += ind + iby*ldha;
    if ( full ) {
        /* full block-column */
#pragma unroll
        for( int j=0; j < BLK_Y; ++j ) {
            tmp = A[j*lda];

            if ( fabs(tmp) > rmax ) {
                cuda_dlag2h_flag = 1;
            }

            HA[j*ldha] = __double2half(tmp);
        }
    }
    else {
        /* partial block-column */
        for( int j=0; (j < BLK_Y) && (iby+j < n); ++j ) {
            tmp = A[j*lda];
            if ( fabs(tmp) > rmax ) {
                cuda_dlag2h_flag = 1;
            }

            HA[j*ldha] = __double2half(tmp);
        }
    }
}


/**
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_dlag2h converts a double-real matrix, A, to a half-real matrix,
 * HA.
 *
 * RMAX is the overflow for the single-real arithmetic.  DLAG2H checks that
 * all the entries of A are between -RMAX and RMAX. If not, the conversion is
 * aborted and a magma_dlag2h_flag is raised.
 *
 *  @param[in] m
 *          The number of lines of the matrix A.  m >= 0.
 *
 *  @param[in] n
 *          The number of columns of the matrix A.  n >= 0.
 *
 *  @param[in] A
 *          On entry, the lda-by-n coefficient matrix A.
 *
 *  @param[in] lda
 *          The leading dimension of the array A.  LDA >= max(1,m).
 *
 *  @param[out] HA
 *          On exit, if INFO=0, the ldha-by-n coefficient matrix HA;
 *          if INFO > 0, the content of HA is unspecified.
 *
 *  @param[in] ldha
 *          The leading dimension of the array HA.  LDHA >= max(1,m).
 *
 *  @param[out]
 *    -     = 0:  successful exit.
 *    -     < 0:  if INFO = -i, the i-th argument had an illegal value
 *    -     = 1:  an entry of the matrix A is greater than the SINGLE PRECISION
 *                overflow threshold, in this case, the content
 *                of HA on exit is unspecified.
 *
 *  @param[in] handle
 *          Cublas handle to execute in.
 *
 **/
extern "C" int
CUDA_dlag2h( int m, int n,
             const double *A, int lda,
             CHAMELEON_Real16_t *HA, int ldha,
             cublasHandle_t handle )
{
    cudaStream_t stream;
    cudaError_t  err;
    double       rmax;

    if ( m < 0 ) {
        return -1;
    }
    else if ( n < 0 ) {
        return -2;
    }
    else if ( lda < chameleon_max(1,m) ) {
        return -4;
    }
    else if ( ldha < chameleon_max(1,m) ) {
        return -6;
    }

    /* quick return */
    if ( m == 0 || n == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    dim3 threads( BLK_X, 1 );
    dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

    /*
     * There is no lapackf77_hlamch, please visit:
     * https://blogs.mathworks.com/cleve/2017/05/08/half-precision-16-bit-floating-point-arithmetic/
     */
    rmax = 65504.;

    cublasGetStream( handle, &stream );

    cuda_dlag2h_kernel<<< grid, threads, 0, stream >>>( m, n, A, lda, HA, ldha, rmax );

    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_dlag2h failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }

    return CHAMELEON_SUCCESS;
}

/*
 * Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
 * Each block has BLK_X threads.
 * Each thread loops across one row, updating BLK_Y entries.
 */
__global__
void cuda_hlag2d_kernel(
    int m, int n,
    const CHAMELEON_Real16_t *HA, int ldha,
    double *A, int lda )
{
    int ind = blockIdx.x * BLK_X + threadIdx.x;
    int iby = blockIdx.y * BLK_Y;

    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);

    /* do only rows inside matrix */
    if ( ind > m ) {
        return;
    }

    A  += ind + iby*lda;
    HA += ind + iby*ldha;
    if ( full ) {
        // full block-column
#pragma unroll
        for( int j=0; j < BLK_Y; ++j ) {
#if defined(PRECISION_zc)
            A[j*ldha] = make_double( HA[j*ldha].x, HA[j*ldha].y );
#else
            A[j*ldha] = HA[j*ldha];
#endif
        }
    }
    else {
        // partial block-column
        for( int j=0; (j < BLK_Y) && (iby+j) < n; ++j ) {
#if defined(PRECISION_zc)
            A[j*ldha] = make_double( HA[j*ldha].x, HA[j*ldha].y );
#else
            A[j*ldha] = HA[j*ldha];
#endif
        }
    }
}


/**
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_hlag2d converts a half-real matrix, HA, to a double-real matrix,
 * A.
 *
 *  Note that while it is possible to overflow while converting from double to
 *  single, it is not possible to overflow when converting from single to
 *  double.
 *
 *  @param[in] m
 *          The number of lines of the matrix A and HA.  m >= 0.
 *
 *  @param[in] n
 *          The number of columns of the matrix A and HA.  n >= 0.
 *
 *  @param[in] HA
 *          On entry, the lda-by-n coefficient matrix HA.
 *
 *  @param[in] ldha
 *          The leading dimension of the array HA.  ldha >= max(1,m).
 *
 *  @param[out] A
 *          On exit, the lda-by-n coefficient matrix A.
 *
 *  @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 *  @param[out]
 *    -     = 0:  successful exit.
 *    -     < 0:  if INFO = -i, the i-th argument had an illegal value
 *
 *  @param[in] handle
 *          Cublas handle to execute in.
 *
 **/
extern "C" int
CUDA_hlag2d( int m, int n,
             const CHAMELEON_Real16_t *HA, int ldha,
             double      *A,  int lda,
             cublasHandle_t handle )
{
    cudaStream_t stream;
    cudaError_t  err;

    if ( m < 0 ) {
        return -1;
    }
    else if ( n < 0 ) {
        return -2;
    }
    else if ( ldha < chameleon_max(1,m) ) {
        return -4;
    }
    else if ( lda < chameleon_max(1,m) ) {
        return -6;
    }

    /* quick return */
    if ( (m == 0) || (n == 0) ) {
        return CHAMELEON_SUCCESS;
    }

    dim3 threads( BLK_X, 1 );
    dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

    cublasGetStream( handle, &stream );
    cuda_hlag2d_kernel<<< grid, threads, 0, stream >>> ( m, n, HA, ldha, A, lda );

    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_hlag2d failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }

    return CHAMELEON_SUCCESS;
}
