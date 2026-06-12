/**
 *
 * @file cuda_zlarfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zlarfb GPU kernel
 *
 * Code originated from MAGMA
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
CUDA_zlarfb( cham_side_t side, cham_trans_t trans,
             cham_dir_t direct, cham_store_t storev,
             int M, int N, int K,
             const cuDoubleComplex *V, int LDV,
             const cuDoubleComplex *T, int LDT,
                   cuDoubleComplex *C, int LDC,
                   cuDoubleComplex *WORK, int LDWORK,
             cublasHandle_t handle )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zzero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
#else
    double zzero = 0.0;
    double zone  = 1.0;
    double mzone = -1.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    cublasStatus_t rc;
    cham_trans_t   transT, notransV, transV;
    cham_uplo_t    uplo;
    int            info = 0;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        return -1;
    }
    if ((trans != ChamNoTrans) && (trans != ChamConjTrans)) {
        return -2;
    }
    if ((direct != ChamDirForward) && (direct != ChamDirBackward)) {
        return -3;
    }
    if ((storev != ChamColumnwise) && (storev != ChamRowwise)) {
        return -4;
    }
    if (M < 0) {
        return -5;
    }
    if (N < 0) {
        return -6;
    }
    if (K < 0) {
        return -9;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0)) {
        return CHAMELEON_SUCCESS;
    }

    // opposite of trans
    if (trans == ChamNoTrans) {
        transT = ChamConjTrans;
    }
    else {
        transT = ChamNoTrans;
    }

    // whether T is upper or lower triangular
    if (direct == ChamDirForward) {
        uplo = ChamUpper;
    }
    else {
        uplo = ChamLower;
    }

    if (storev == ChamColumnwise) {
        notransV = ChamNoTrans;
        transV   = ChamConjTrans;
    }
    else {
        notransV = ChamConjTrans;
        transV   = ChamNoTrans;
    }

    if ( side == ChamLeft ) {
        // Form H C or H^H C
        // Comments assume H C. When forming H^H C, T gets transposed via transT.

        // W = C^H V
        rc = cublasZgemm( handle,
                          chameleon_cublas_const(ChamConjTrans), chameleon_cublas_const(notransV),
                          N, K, M,
                          CUBLAS_SADDR(zone),  C, LDC,
                                               V, LDV,
                          CUBLAS_SADDR(zzero), WORK, LDWORK );
        assert( rc == CUBLAS_STATUS_SUCCESS );
        info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

        // W = W T^H = C^H V T^H
        CUDA_ztrmm( ChamRight, uplo, transT, ChamNonUnit,
                    N, K,
                    &zone, T,    LDT,
                           WORK, LDWORK,
                    handle );

        // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
        rc = cublasZgemm( handle,
                          chameleon_cublas_const(notransV), chameleon_cublas_const(ChamConjTrans),
                          M, N, K,
                          CUBLAS_SADDR(mzone), V,    LDV,
                                               WORK, LDWORK,
                          CUBLAS_SADDR(zone),  C,    LDC );
        assert( rc == CUBLAS_STATUS_SUCCESS );
        info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
    }
    else {
        // Form C H or C H^H
        // Comments assume C H. When forming C H^H, T gets transposed via trans.

        // W = C V
        rc = cublasZgemm( handle,
                          chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(notransV),
                          M, K, N,
                          CUBLAS_SADDR(zone),  C, LDC,
                                               V, LDV,
                          CUBLAS_SADDR(zzero), WORK, LDWORK );
        assert( rc == CUBLAS_STATUS_SUCCESS );
        info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

        // W = W T = C V T
        CUDA_ztrmm( ChamRight, uplo, trans, ChamNonUnit,
                    M, K,
                    &zone, T,    LDT,
                           WORK, LDWORK,
                    handle );

        // C = C - W V^H = C - C V T V^H = C (I - V T V^H) = C H
        rc = cublasZgemm( handle,
                          chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(transV),
                          M, N, K,
                          CUBLAS_SADDR(mzone), WORK, LDWORK,
                                               V,    LDV,
                          CUBLAS_SADDR(zone),  C,    LDC );
        assert( rc == CUBLAS_STATUS_SUCCESS );
        info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
    }

    return (info == 0) ? CHAMELEON_SUCCESS : CHAMELEON_ERR_UNEXPECTED;
}
