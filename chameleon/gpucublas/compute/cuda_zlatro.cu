/**
 *
 * @file cuda_zlatro.cu
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlatro GPU kernel
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2026-01-23
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#if defined( PRECISION_z )
#define cudaConj( val ) cuConj( val )
#elif defined( PRECISION_c )
#define cudaConj( val ) cuConjf( val )
#else
#define cudaConj( val ) val
#endif

#define BLK_X 32
#define BLK_Y 32

__global__
void CUDA_zlatro_upper_trans( int m, int n,
                              const cuDoubleComplex *A, int lda,
                              cuDoubleComplex *B, int ldb )
{
    __shared__ cuDoubleComplex wA[BLK_X * BLK_Y];

    int i = threadIdx.x;
    int j = threadIdx.y;

    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /* In Lower part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if( ( x > y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    wA[ j + BLK_Y * i ] = A[ x + y*lda ];
    __syncthreads();
    B[ y + ldb * x ] = wA[ j + BLK_Y * i ];
}

__global__
void CUDA_zlatro_lower_trans( int m, int n,
                              const cuDoubleComplex *A, int lda,
                              cuDoubleComplex *B, int ldb )
{
    __shared__ cuDoubleComplex wA[BLK_X * BLK_Y];

    int i = threadIdx.x;
    int j = threadIdx.y;

    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /* In Upper part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if( ( x < y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    wA[ j + BLK_Y * i ] = A[ x + y*lda ];
    __syncthreads();
    B[ y + ldb * x ] = wA[ j + BLK_Y * i ];
}

#if defined(PRECISION_z) || defined(PRECISION_c)
__global__
void CUDA_zlatro_upper_conjtrans( int m, int n,
                                  const cuDoubleComplex *A, int lda,
                                  cuDoubleComplex *B, int ldb )
{
    __shared__ cuDoubleComplex wA[BLK_X * BLK_Y];

    int i = threadIdx.x;
    int j = threadIdx.y;

    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /* In Lower part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if( ( x > y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    wA[ j + BLK_Y * i ] = A[ x + y*lda ];
    __syncthreads();
    B[ y + ldb * x ] = cudaConj(wA[ j + BLK_Y * i ]);
}

__global__
void CUDA_zlatro_lower_conjtrans( int m, int n,
                                  const cuDoubleComplex *A, int lda,
                                  cuDoubleComplex *B, int ldb )
{
    __shared__ cuDoubleComplex wA[BLK_X * BLK_Y];

    int i = threadIdx.x;
    int j = threadIdx.y;

    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /* In Upper part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if( ( x < y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    wA[ j + BLK_Y * i ] = A[ x + y*lda ];
    __syncthreads();
    B[ y + ldb * x ] = cudaConj(wA[ j + BLK_Y * i ]);
}
#endif

extern "C" int
CUDA_zlatro( cham_uplo_t uplo, cham_trans_t trans,
             int m, int n,
             const cuDoubleComplex *A, int lda,
             cuDoubleComplex *B, int ldb,
             cublasHandle_t handle )
{
    cudaError_t err;
    cublasStatus_t rc = CUBLAS_STATUS_SUCCESS;
    cudaStream_t stream;

#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex one = make_cuDoubleComplex(1., 0.);
    cuDoubleComplex zero = make_cuDoubleComplex(0., 0.);
#else
    double one = 1.;
    double zero = 0.;
#endif

    dim3 threads(BLK_X, BLK_Y); /* Distribution by column */
    dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

    cublasGetStream( handle, &stream );

    switch ( uplo ) {
        case ChamUpperLower:
            rc = cublasZgeam( handle,
                              (cublasOperation_t)chameleon_cublas_const(trans),
                              (cublasOperation_t)chameleon_cublas_const(ChamNoTrans),
                              n, m,
                              &one, A, lda,
                              &zero, B, ldb,
                              B, ldb );
            break;
        case ChamUpper:
#if defined(PRECISION_z) || defined(PRECISION_c)
            if( trans == ChamConjTrans ) {
                CUDA_zlatro_upper_conjtrans<<<grid, threads, 0, stream>>>( m, n,
                                                                           A, lda,
                                                                           B, ldb );
            }
            else
#endif
            {
                CUDA_zlatro_upper_trans<<<grid, threads, 0, stream>>>( m, n,
                                                                       A, lda,
                                                                       B, ldb );
            }
            break;
        case ChamLower:
#if defined(PRECISION_z) || defined(PRECISION_c)
            if( trans == ChamConjTrans ) {
                CUDA_zlatro_lower_conjtrans<<<grid, threads, 0, stream>>>( m, n,
                                                                           A, lda,
                                                                           B, ldb );
            }
            else
#endif
            {
                CUDA_zlatro_lower_trans<<<grid, threads, 0, stream>>>( m, n,
                                                                       A, lda,
                                                                       B, ldb );
            }
            break;
    }

    assert( rc == CUBLAS_STATUS_SUCCESS );
    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_zlatro failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }
    return rc;
}
