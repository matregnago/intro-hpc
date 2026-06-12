/**
 *
 * @file cuda_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlacpy GPU kernel
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#define BLK_X 32
#define BLK_Y 32

__global__
void CUDA_zlacpy_upper( int m, int n,
                        const cuDoubleComplex *A, int lda,
                        cuDoubleComplex *B, int ldb )
{
    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /*
     * In Lower part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if ( ( x > y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    B[ ldb * y + x ] = A[ lda * y + x ];
}

__global__
void CUDA_zlacpy_lower( int m, int n,
                        const cuDoubleComplex *A, int lda,
                        cuDoubleComplex *B, int ldb )
{
    int x = threadIdx.x + blockIdx.x * BLK_X;
    int y = threadIdx.y + blockIdx.y * BLK_Y;

    /*
     * In Upper part of A - quick return
     * We are launching lots of blocks for nothing :/
     */
    if ( ( x < y ) || ( x >= m ) || ( y >= n ) ) {
        return;
    }

    B[ ldb * y + x ] = A[ lda * y + x ];
}

extern "C" int
CUDA_zlacpy( cham_uplo_t            uplo,
             int                    m,
             int                    n,
             const cuDoubleComplex *A,
             int                    lda,
             cuDoubleComplex       *B,
             int                    ldb,
             cublasHandle_t         handle )
{
    cudaError_t    err;
    cublasStatus_t rc = CUBLAS_STATUS_SUCCESS;
    cudaStream_t   stream;

#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex one  = make_cuDoubleComplex(1., 0.);
    cuDoubleComplex zero = make_cuDoubleComplex(0., 0.);
#else
    double one  = 1.;
    double zero = 0.;
#endif

    dim3 threads( BLK_X, BLK_Y );
    dim3 grid( chameleon_ceil( m, BLK_X ), chameleon_ceil( n, BLK_Y ) );

    cublasGetStream( handle, &stream );

    switch ( uplo ) {
        case ChamUpperLower:
            /* Use GEAM as recommended by CUDA API */
            rc = cublasZgeam( handle,
                              (cublasOperation_t)chameleon_cublas_const(ChamNoTrans),
                              (cublasOperation_t)chameleon_cublas_const(ChamNoTrans),
                              m, n,
                              &one,  A, lda,
                              &zero, B, ldb,
                              B, ldb );
            break;
        case ChamUpper:
            {
                CUDA_zlacpy_upper<<<grid, threads, 0, stream>>>( m, n,
                                                                 A, lda,
                                                                 B, ldb );
            }
            break;
        case ChamLower:
            {
                CUDA_zlacpy_lower<<<grid, threads, 0, stream>>>( m, n,
                                                                 A, lda,
                                                                 B, ldb );
            }
            break;
    }

    assert( rc == CUBLAS_STATUS_SUCCESS );
    err = cudaGetLastError();
    if ( err != cudaSuccess )
    {
        fprintf( stderr, "CUDA_zlacpy failed to launch CUDA kernel %s\n", cudaGetErrorString(err) );
        return CHAMELEON_ERR_UNEXPECTED;
    }
    return rc;
}
