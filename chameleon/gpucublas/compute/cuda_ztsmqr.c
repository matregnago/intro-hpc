/**
 *
 * @file cuda_ztsmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztsmqr GPU kernel
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
CUDA_ztsmqr( cham_side_t side, cham_trans_t trans,
             int M1, int N1,
             int M2, int N2,
             int K, int IB,
             cuDoubleComplex *A1,    int LDA1,
             cuDoubleComplex *A2,    int LDA2,
             const cuDoubleComplex *V,     int LDV,
             const cuDoubleComplex *T,     int LDT,
             cuDoubleComplex *WORK,  int LWORK,
             cublasHandle_t handle )
{
    int i, i1, i3;
    int NQ;
    int kb;
    int ic = 0;
    int jc = 0;
    int mi = M1;
    int ni = N1;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        return -1;
    }

    /* NQ is the order of Q */
    if (side == ChamLeft) {
        NQ = M2;
    }
    else {
        NQ = N2;
    }

    if ((trans != ChamNoTrans) && (trans != ChamConjTrans)) {
        return -2;
    }
    if (M1 < 0) {
        return -3;
    }
    if (N1 < 0) {
        return -4;
    }
    if ( (M2 < 0) ||
         ( (M2 != M1) && (side == ChamRight) ) ){
        return -5;
    }
    if ( (N2 < 0) ||
         ( (N2 != N1) && (side == ChamLeft) ) ){
        return -6;
    }
    if ((K < 0) ||
        ( (side == ChamLeft)  && (K > M1) ) ||
        ( (side == ChamRight) && (K > N1) ) ) {
        return -7;
    }
    if (IB < 0) {
        return -8;
    }
    if (LDA1 < chameleon_max(1,M1)){
        return -10;
    }
    if (LDA2 < chameleon_max(1,M2)){
        return -12;
    }
    if (LDV < chameleon_max(1,NQ)){
        return -14;
    }
    if (LDT < chameleon_max(1,IB)){
        return -16;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( ((side == ChamLeft ) && (trans != ChamNoTrans))  ||
         ((side == ChamRight) && (trans == ChamNoTrans)) )
    {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    for (i = i1; (i > -1) && (i < K); i+=i3) {
        kb = chameleon_min(IB, K-i);

        if (side == ChamLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M1 - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N1 - i;
            jc = i;
        }

        /*
         * Apply H or H' (NOTE: CORE_zparfb used to be CORE_ztsrfb)
         */
        CUDA_zparfb(
            side, trans, ChamDirForward, ChamColumnwise,
            mi, ni, M2, N2, kb, 0,
            A1 + LDA1*jc+ic, LDA1,
            A2, LDA2,
            V + LDV*i, LDV,
            T + LDT*i, LDT,
            WORK, LWORK, handle );
    }
    return CHAMELEON_SUCCESS;
}
