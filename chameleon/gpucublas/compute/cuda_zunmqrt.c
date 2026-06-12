/**
 *
 * @file cuda_zunmqrt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zunmqrt GPU kernel
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
CUDA_zunmqrt( cham_side_t side, cham_trans_t trans,
              int M, int N, int K, int IB,
              const cuDoubleComplex *A,    int LDA,
              const cuDoubleComplex *T,    int LDT,
              cuDoubleComplex *C,    int LDC,
              cuDoubleComplex *WORK, int LDWORK,
              cublasHandle_t handle )
{
    int i, kb;
    int i1, i3;
    int nq, nw;
    int ic = 0;
    int jc = 0;
    int ni = N;
    int mi = M;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        return -1;
    }
    /*
     * NQ is the order of Q and NW is the minimum dimension of WORK
     */
    if (side == ChamLeft) {
        nq = M;
        nw = N;
    }
    else {
        nq = N;
        nw = M;
    }

    if ((trans != ChamNoTrans) && (trans != ChamConjTrans)) {
        return -2;
    }
    if (M < 0) {
        return -3;
    }
    if (N < 0) {
        return -4;
    }
    if ((K < 0) || (K > nq)) {
        return -5;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        return -6;
    }
    if ((LDA < chameleon_max(1,nq)) && (nq > 0)) {
        return -8;
    }
    if ((LDC < chameleon_max(1,M)) && (M > 0)) {
        return -12;
    }
    if ((LDWORK < chameleon_max(1,nw)) && (nw > 0)) {
        return -14;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0))
        return CHAMELEON_SUCCESS;

    if (((side == ChamLeft) && (trans != ChamNoTrans))
        || ((side == ChamRight) && (trans == ChamNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    for(i = i1; (i >- 1) && (i < K); i+=i3 ) {
        kb = chameleon_min(IB, K-i);

        if (side == ChamLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N - i;
            jc = i;
        }

        CUDA_zlarfb( side, trans, ChamDirForward, ChamColumnwise,
                     mi, ni, kb,
                     A + LDA * i  + i,  LDA,
                     T + LDT * i,       LDT,
                     C + LDC * jc + ic, LDC,
                     WORK, LDWORK,
                     handle );
    }

    return CHAMELEON_SUCCESS;
}
