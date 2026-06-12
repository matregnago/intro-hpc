/**
 *
 * @file core_zttmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zttmlq CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zttmlq overwrites the general complex M1-by-N1 tile A1 and
 *  M2-by-N2 tile A2 (N1 == N2) with
 *
 *                        SIDE = 'L'        SIDE = 'R'
 *    TRANS = 'N':         Q * | A1 |       | A1 | * Q
 *                             | A2 |       | A2 |
 *
 *    TRANS = 'C':      Q^H * | A1 |       | A1 | * Q^H
 *                             | A2 |       | A2 |
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zttqrt.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q^H from the Left;
 *         @arg ChamRight : apply Q or Q^H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  ConjTranspose, apply Q^H.
 *
 * @param[in] M1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2. M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTTQRT in the first k rows of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[out] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size LDWORK-by-N1.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK. LDWORK >= max(1,IB).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zttmlq(cham_side_t side, cham_trans_t trans,
                int M1, int N1, int M2, int N2, int K, int IB,
                CHAMELEON_Complex64_t *A1, int LDA1,
                CHAMELEON_Complex64_t *A2, int LDA2,
                const CHAMELEON_Complex64_t *V, int LDV,
                const CHAMELEON_Complex64_t *T, int LDT,
                CHAMELEON_Complex64_t *WORK, int LDWORK)
{
    int i, i1, i3;
    int NW;
    int kb, l;
    int ic = 0;
    int jc = 0;
    int mi1 = M1;
    int mi2 = M2;
    int ni1 = N1;
    int ni2 = N2;

    /* Check input arguments */
    if ((side != ChamLeft) && (side != ChamRight)) {
        coreblas_error(1, "Illegal value of side");
        return -1;
    }

    /* NQ is the order of Q */
    if (side == ChamLeft) {
        NW = IB;
    }
    else {
        NW = N1;
    }

    if ((trans != ChamNoTrans) && (trans != ChamConjTrans)) {
        coreblas_error(2, "Illegal value of trans");
        return -2;
    }
    if (M1 < 0) {
        coreblas_error(3, "Illegal value of M1");
        return -3;
    }
    if (N1 < 0) {
        coreblas_error(4, "Illegal value of N1");
        return -4;
    }
    if ((M2 < 0) ||
        ( (side == ChamRight) && (M1 != M2) ) ) {
        coreblas_error(5, "Illegal value of M2");
        return -5;
    }
    if ((N2 < 0) ||
        ( (side == ChamLeft) && (N1 != N2) ) ) {
        coreblas_error(6, "Illegal value of N2");
        return -6;
    }
    if ((K < 0) ||
        ( (side == ChamLeft)  && (K > M1) ) ||
        ( (side == ChamRight) && (K > N1) ) ) {
        coreblas_error(7, "Illegal value of K");
        return -7;
    }
    if (IB < 0) {
        coreblas_error(8, "Illegal value of IB");
        return -8;
    }
    if (LDA1 < chameleon_max(1,M1)){
        coreblas_error(10, "Illegal value of LDA1");
        return -10;
    }
    if (LDA2 < chameleon_max(1,M2)){
        coreblas_error(12, "Illegal value of LDA2");
        return -12;
    }
    if (LDV < chameleon_max(1,K)){
        coreblas_error(14, "Illegal value of LDV");
        return -14;
    }
    if (LDT < chameleon_max(1,IB)){
        coreblas_error(16, "Illegal value of LDT");
        return -16;
    }
    if (LDWORK < chameleon_max(1,NW)){
        coreblas_error(18, "Illegal value of LDWORK");
        return -18;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( ((side == ChamLeft ) && (trans == ChamNoTrans)) ||
         ((side == ChamRight) && (trans != ChamNoTrans)) )
    {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    /* Transpose */
    if (trans == ChamNoTrans) {
        trans = ChamConjTrans;
    }
    else {
        trans = ChamNoTrans;
    }

    for (i = i1; (i > -1) && (i < K); i+=i3) {
        kb = chameleon_min(IB, K-i);

        if (side == ChamLeft) {
            mi1 = kb; // M1 - i;
            mi2 = chameleon_min(i+kb, M2);
            l   = chameleon_min(kb, chameleon_max(0, M2-i));
            ic  = i;
        }
        else {
            ni1 = kb;
            ni2 = chameleon_min(i+kb, N2);
            l   = chameleon_min(kb, chameleon_max(0, N2-i));
            jc  = i;
        }

        /*
         * Apply H or H'
         */
        CORE_zparfb(
            side, trans, ChamDirForward, ChamRowwise,
            mi1, ni1, mi2, ni2, kb, l,
            A1 + LDA1 * jc + ic, LDA1,
            A2, LDA2,
            V + i,       LDV,
            T + LDT * i, LDT,
            WORK, LDWORK );
    }
    return CHAMELEON_SUCCESS;
}
