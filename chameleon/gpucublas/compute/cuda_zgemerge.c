/**
 *
 * @file cuda_zgemerge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgemerge GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

int
CUDA_zgemerge( cham_side_t side, cham_diag_t diag,
               int M, int N,
               const cuDoubleComplex *A, int LDA,
               cuDoubleComplex *B, int LDB,
               cublasHandle_t handle )
{
    const cuDoubleComplex *cola;
    cuDoubleComplex       *colb;
    cublasStatus_t         rc;
    int                    i;

    if (M < 0) {
        return -1;
    }
    if (N < 0) {
        return -2;
    }
    if ( (LDA < chameleon_max(1,M)) && (M > 0) ) {
        return -5;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        return -7;
    }

    if (side == ChamLeft){
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
            rc = cublasZcopy( handle, i+1, cola, 1, colb, 1 );
            assert( rc == CUBLAS_STATUS_SUCCESS );
        }
    }else{
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
            rc = cublasZcopy( handle, M-i, cola, 1, colb, 1 );
            assert( rc == CUBLAS_STATUS_SUCCESS );
        }
    }

    (void)diag;
    (void)rc;
    return CHAMELEON_SUCCESS;
}
