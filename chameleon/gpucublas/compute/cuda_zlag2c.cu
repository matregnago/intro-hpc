/**
 *
 * @file cuda_zlag2c.cu
 *
 * @copyright 2023-2023 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlag2c GPU kernel
 *
 * @version 1.4.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions mixed zc -> ds
 *
 * This file is an adaptation of the MAGMA zlag2c.cu and clag2z files.
 *
 */
#include "gpucublas.h"
#include <coreblas/lapacke.h>

#define BLK_X 64
#define BLK_Y 32

__device__  int cuda_zlag2c_flag = 0;

/*
 * Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
 * Each block has BLK_X threads.
 * Each thread loops across one row, updating BLK_Y entries.
 */
__global__
void cuda_zlag2c_kernel(
    int m, int n,
    const cuDoubleComplex *A,  int lda,
    cuFloatComplex        *SA, int ldsa,
    double rmax )
{
    cuDoubleComplex tmp;

    int ind = blockIdx.x * BLK_X + threadIdx.x;
    int iby = blockIdx.y * BLK_Y;

    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);

    /* do only rows inside matrix */
    if ( ind > m ) {
        return;
    }

    A  += ind + iby*lda;
    SA += ind + iby*ldsa;
    if ( full ) {
        /* full block-column */
#pragma unroll
        for( int j=0; j < BLK_Y; ++j ) {
            tmp = A[j*lda];

            if (
#if defined(PRECISION_zc)
                (fabs(tmp.x) > rmax) || (fabs(tmp.y) > rmax)
#else
                (fabs(tmp) > rmax)
#endif
                )
            {
                cuda_zlag2c_flag = 1;
            }

#if defined(PRECISION_zc)
            SA[j*ldsa] = make_cuFloatComplex(tmp.x, tmp.y);
#else
            SA[j*ldsa] = tmp;
#endif
        }
    }
    else {
        /* partial block-column */
        for( int j=0; (j < BLK_Y) && (iby+j < n); ++j ) {
            tmp = A[j*lda];
            if (
#if defined(PRECISION_zc)
                (fabs(tmp.x) > rmax) || (fabs(tmp.y) > rmax)
#else
                (fabs(tmp) > rmax)
#endif
                )
            {
                cuda_zlag2c_flag = 1;
            }

#if defined(PRECISION_zc)
            SA[j*ldsa] = make_cuFloatComplex(tmp.x, tmp.y);
#else
            SA[j*ldsa] = tmp;
#endif
        }
    }
}


/**
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_zlag2c converts a double-complex matrix, A, to a single-complex matrix,
 * SA.
 *
 * RMAX is the overflow for the single-complex arithmetic.  ZLAG2C checks that
 * all the entries of A are between -RMAX and RMAX. If not, the conversion is
 * aborted and a magma_zlag2c_flag is raised.
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
 *  @param[out] SA
 *          On exit, if INFO=0, the ldsa-by-n coefficient matrix SA;
 *          if INFO > 0, the content of SA is unspecified.
 *
 *  @param[in] ldsa
 *          The leading dimension of the array SA.  LDSA >= max(1,m).
 *
 *  @param[out]
 *    -     = 0:  successful exit.
 *    -     < 0:  if INFO = -i, the i-th argument had an illegal value
 *    -     = 1:  an entry of the matrix A is greater than the COMPLEX
 *                overflow threshold, in this case, the content
 *                of SA on exit is unspecified.
 *
 *  @param[in] handle
 *          Cublas handle to execute in.
 *
 **/
extern "C" int
CUDA_zlag2c( int m, int n,
             const cuDoubleComplex *A, int lda,
             cuFloatComplex *SA, int ldsa,
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
    else if ( ldsa < chameleon_max(1,m) ) {
        return -6;
    }

    /* quick return */
    if ( m == 0 || n == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    dim3 threads( BLK_X, 1 );
    dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

    rmax = LAPACKE_slamch_work( 'O' );
    cublasGetStream( handle, &stream );

    cuda_zlag2c_kernel<<< grid, threads, 0, stream >>>( m, n, A, lda, SA, ldsa, rmax );

    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_zlag2c failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
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
void cuda_clag2z_kernel(
    int m, int n,
    const cuFloatComplex *SA, int ldsa,
    cuDoubleComplex       *A, int lda )
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
    SA += ind + iby*ldsa;
    if ( full ) {
        // full block-column
#pragma unroll
        for( int j=0; j < BLK_Y; ++j ) {
#if defined(PRECISION_zc)
            A[j*ldsa] = make_cuDoubleComplex( SA[j*ldsa].x, SA[j*ldsa].y );
#else
            A[j*ldsa] = SA[j*ldsa];
#endif
        }
    }
    else {
        // partial block-column
        for( int j=0; (j < BLK_Y) && (iby+j) < n; ++j ) {
#if defined(PRECISION_zc)
            A[j*ldsa] = make_cuDoubleComplex( SA[j*ldsa].x, SA[j*ldsa].y );
#else
            A[j*ldsa] = SA[j*ldsa];
#endif
        }
    }
}


/**
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_clag2z converts a single-complex matrix, SA, to a double-complex matrix,
 * A.
 *
 *  Note that while it is possible to overflow while converting from double to
 *  single, it is not possible to overflow when converting from single to
 *  double.
 *
 *  @param[in] m
 *          The number of lines of the matrix A and SA.  m >= 0.
 *
 *  @param[in] n
 *          The number of columns of the matrix A and SA.  n >= 0.
 *
 *  @param[in] SA
 *          On entry, the lda-by-n coefficient matrix SA.
 *
 *  @param[in] ldsa
 *          The leading dimension of the array SA.  ldsa >= max(1,m).
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
CUDA_clag2z( int m, int n,
             const cuFloatComplex *SA, int ldsa,
             cuDoubleComplex      *A,  int lda,
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
    else if ( ldsa < chameleon_max(1,m) ) {
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
    cuda_clag2z_kernel<<< grid, threads, 0, stream >>> ( m, n, SA, ldsa, A, lda );

    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_clag2z failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }

    return CHAMELEON_SUCCESS;
}
