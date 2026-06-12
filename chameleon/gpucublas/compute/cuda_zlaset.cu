/**
 *
 * @file cuda_zlaset.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800 ), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlaset GPU kernel
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2026-02-24
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#define BLK_X 32
#define BLK_Y 32

__global__
void CUDA_zlaset_upper( int              m,
                        int              n,
                        cuDoubleComplex  alpha,
                        cuDoubleComplex  beta,
                        cuDoubleComplex *A,
                        int              lda )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if ( ( i >= m ) || ( j >= n ) || ( j < i ) ) {
        return;
    }

    if ( j == i ) {
        A[i + j*lda] = beta;
    }
    else {
        A[i + j*lda] = alpha;
    }
}

__global__
void CUDA_zlaset_lower( int              m,
                        int              n,
                        cuDoubleComplex  alpha,
                        cuDoubleComplex  beta,
                        cuDoubleComplex *A,
                        int              lda )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if ( ( i >= m ) || ( j >= n ) || ( i < j ) ) {
        return;
    }

    if ( j == i ) {
        A[i + j*lda] = beta;
    }
    else {
        A[i + j*lda] = alpha;
    }
}

__global__
void CUDA_zlaset_upperlower( int              m,
                             int              n,
                             cuDoubleComplex  alpha,
                             cuDoubleComplex  beta,
                             cuDoubleComplex *A,
                             int              lda )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if ( ( i >= m ) || ( j >= n ) ) {
        return;
    }
    if ( j == i ) {
        A[i + j*lda] = beta;
    }
    else {
        A[i + j*lda] = alpha;
    }
}

extern "C" int
CUDA_zlaset( cham_uplo_t uplo, int m, int n,
             const cuDoubleComplex *alpha, const cuDoubleComplex *beta,
             cuDoubleComplex *A, int lda,
             cublasHandle_t handle )
{
    cudaError_t err;
    cudaStream_t stream;
    struct cudaPointerAttributes attr;

    cublasGetStream( handle, &stream );

    cudaPointerGetAttributes( &attr, alpha );
    assert( attr.memoryType == cudaMemoryTypeHost );
    assert( attr.hostPointer == alpha );
    assert( attr.devicePointer == NULL );

    cudaPointerGetAttributes( &attr, beta );
    assert( attr.memoryType == cudaMemoryTypeHost );
    assert( attr.hostPointer == beta );
    assert( attr.devicePointer == NULL );

#if defined(PRECISION_z) || defined(PRECISION_c)
    if( ( cuCreal(*alpha) == 0. ) && ( cuCimag(*alpha) == 0. ) &&
        ( cuCreal(*beta)  == 0. ) && ( cuCimag(*beta)  == 0. ) &&
        ( uplo == ChamUpperLower ) )
#else
    if( ( *alpha == 0. ) && ( *beta == 0. ) && ( uplo == ChamUpperLower ) )
#endif
    {
        err = cudaMemset2DAsync( A, sizeof(cuDoubleComplex) * lda, 0.,
                                 sizeof(cuDoubleComplex) * m, n, stream );
    }
    else
    {
        dim3 threads(BLK_X, BLK_Y);
        dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

        cublasGetStream( handle, &stream );

        switch ( uplo ) {
            case ChamUpperLower:
                CUDA_zlaset_upperlower<<<grid, threads, 0, stream>>>( m, n,
                                                                      *alpha, *beta,
                                                                      A, lda );
                break;
            case ChamUpper:
                CUDA_zlaset_upper<<<grid, threads, 0, stream>>>( m, n,
                                                                 *alpha, *beta,
                                                                 A, lda );
                break;
            case ChamLower:
                CUDA_zlaset_lower<<<grid, threads, 0, stream>>>( m, n,
                                                                 *alpha, *beta,
                                                                 A, lda );
                break;
        }
    }


    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_zlaset failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }
    return CHAMELEON_SUCCESS;
}
