/**
 *
 * @file cuda_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherfb GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

int
CUDA_zherfb( cham_uplo_t uplo, int n,
             int k, int ib, int nb,
             const cuDoubleComplex *A, int lda,
             const cuDoubleComplex *T, int ldt,
             cuDoubleComplex *C, int ldc,
             cuDoubleComplex *WORK, int ldwork,
             cublasHandle_t handle )
{
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        gpucublas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        gpucublas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        gpucublas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        gpucublas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        gpucublas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < chameleon_max(1,n)) && (n > 0) ) {
        gpucublas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < chameleon_max(1,ib)) && (ib > 0) ) {
        gpucublas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < chameleon_max(1,n)) && (n > 0) ) {
        gpucublas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == ChamLower) {
        /* Left */
        CUDA_zunmqrt( ChamLeft, ChamConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      handle );
        /* Right */
        CUDA_zunmqrt( ChamRight, ChamNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      handle );
    }
    else {
        /* Right */
        CUDA_zunmlqt( ChamRight, ChamConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      handle );
        /* Left */
        CUDA_zunmlqt( ChamLeft, ChamNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      handle );
    }
    return CHAMELEON_SUCCESS;
}
