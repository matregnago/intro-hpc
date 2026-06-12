/**
 *
 * @file cuda_zparfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zparfb GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-01-29
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

/**
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 *  CUDA_zparfb applies a complex upper triangular block reflector H
 *  or its transpose H' to a complex rectangular matrix formed by
 *  coupling two tiles A1 and A2. Matrix V is:
 *
 *          COLUMNWISE                    ROWWISE
 *
 *         |     K     |                 |      N2-L     |   L  |
 *      __ _____________ __           __ _________________        __
 *         |    |      |                 |               | \
 *         |    |      |                 |               |   \    L
 *    M2-L |    |      |              K  |_______________|_____\  __
 *         |    |      | M2              |                      |
 *      __ |____|      |                 |                      | K-L
 *         \    |      |              __ |______________________| __
 *       L   \  |      |
 *      __     \|______| __              |          N2          |
 *
 *         | L |  K-L  |
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q^H from the Left;
 *         @arg ChamRight : apply Q or Q^H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   : No transpose, apply Q;
 *         @arg ChamConjTrans : ConjTranspose, apply Q^H.
 *
 * @param[in] direct
 *         Indicates how H is formed from a product of elementary
 *         reflectors
 *         @arg ChamDirForward  : H = H(1) H(2) . . . H(k) (Forward)
 *         @arg ChamDirBackward : H = H(k) . . . H(2) H(1) (Backward)
 *
 * @param[in] storev
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored:
 *         @arg ChamColumnwise
 *         @arg ChamRowwise
 *
 * @param[in] M1
 *         The number of columns of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of rows of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of columns of the tile A2. M2 >= 0.
 *
 * @param[in] N2
 *         The number of rows of the tile A2. N2 >= 0.
 *
 * @param[in] K
 *         The order of the matrix T (= the number of elementary
 *         reflectors whose product defines the block reflector).
 *
 * @param[in] L
 *         The size of the triangular part of V
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,N2).
 *
 * @param[in] V
 *         (LDV,K) if STOREV = 'C'
 *         (LDV,M2) if STOREV = 'R' and SIDE = 'L'
 *         (LDV,N2) if STOREV = 'R' and SIDE = 'R'
 *         Matrix V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *         If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M2);
 *         if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N2);
 *         if STOREV = 'R', LDV >= K.
 *
 * @param[out] T
 *         The triangular K-by-K matrix T in the representation of the
 *         block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= K.
 *
 * @param[in,out] WORK
 *         Workspace of dimension at least:
 *            - K * (M2 + N2).
 *         If L > 0, it is recommended to extend it to
 *            - K * (2 * M2 + N2 ) if side == ChamLeft.
 *            - K * (M2 + 2 * N2 ) if side == ChamRight.
 *
 * @param[in] LWORK
 *         The dimension of the array WORK. If LWORK < 0, returns immediately
 *         the recommended workspace size.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval  <0 if -i, the i-th argument had an illegal value
 * @retval  The recommended LWORK value, if LWORK == -1 on entry.
 */
/* This kernel is never traced so return type on previous line for convert2eztrace.pl script */
int
CUDA_zparfb( cham_side_t side, cham_trans_t trans,
             cham_dir_t direct, cham_store_t storev,
             int M1, int N1, int M2, int N2, int K, int L,
             cuDoubleComplex       *A1,   int LDA1,
             cuDoubleComplex       *A2,   int LDA2,
             const cuDoubleComplex *V,    int LDV,
             const cuDoubleComplex *T,    int LDT,
             cuDoubleComplex       *WORK, int LWORK,
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

    cuDoubleComplex *workW, *workC, *workV;
    cublasStatus_t   rc;
    cudaStream_t     stream;
    int              i, j, ldW, ldC, ldV;
    cham_trans_t     transW, transA2;
    int              wssize = 0;
    int              wrsize = 0;
    int              info = 0;

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
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == ChamRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == ChamLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }

    if (direct == ChamDirForward) {
        wssize = K * (M2 + N2);
        wrsize = wssize;
        if ( L > 0 ) {
            wrsize +=  (side == ChamLeft) ? M2 * K : K * N2;
        }
    }

    if ( LWORK < 0 ) {
        return wrsize;
    }
    else if ( LWORK < wssize ) {
        gpucublas_error(20, "Illegal value of LWORK");
        return -20;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0)) {
        return CHAMELEON_SUCCESS;
    }

    cublasGetStream( handle, &stream );

    if (direct == ChamDirForward) {

        if (side == ChamLeft) {
            /*
             * Column or Rowwise / Forward / Left
             * ----------------------------------
             *
             * Form  H * A  or  H' * A  where  A = ( A1 )
             *                                     ( A2 )
             */

            /*
             * Store in WORK (N1 == N2):
             *    - Workspace W for the copy of A1 + V' * A2 (K  x N1)
             *    - Workspace C for the copy of V * T        (M2 x K )
             *    - Workspace V for the copy of V            (M2 x K )
             */
            workW = WORK;
            ldW = K;

            workC = workW + K * N1;
            ldC = M2;

            if ( L == 0 ) {
                workV = (cuDoubleComplex*)V;
                ldV   = LDV;
            }
            else {
                if ( LWORK < wrsize ) {
                    workC = NULL;
                    workV = workW + K * N1;
                }
                else {
                    workV = workC + M2 * K;
                }

                if ( storev == ChamColumnwise ) {
                    ldV = M2;

                    /*
                     * Backup V, and put 0 in the lower part
                     */
                    cudaMemcpy2DAsync( workV, sizeof(cuDoubleComplex) * ldV,
                                       V,     sizeof(cuDoubleComplex) * LDV,
                                       sizeof(cuDoubleComplex) * M2, K,
                                       cudaMemcpyDeviceToDevice, stream );

                    j = 0;
                    i = M2 - L + 1;
                    for(; (i < M2) && (j < K); i++, j++ ) {
                        cudaMemsetAsync( workV + j * ldV + i, 0,
                                         sizeof(cuDoubleComplex) * (M2 - i),
                                         stream );
                    }
                }
                else {
                    ldV = K;

                    /*
                     * Backup V, and put 0 in the upper part
                     */
                    cudaMemcpy2DAsync( workV, sizeof(cuDoubleComplex) * ldV,
                                       V,     sizeof(cuDoubleComplex) * LDV,
                                       sizeof(cuDoubleComplex) * K, M2,
                                       cudaMemcpyDeviceToDevice, stream );

                    for(j = 1; j < K; j++) {
                        cudaMemsetAsync( workV + ldV * ( M2 - L + j ), 0,
                                         sizeof(cuDoubleComplex) * j,
                                         stream );
                    }
                }
            }

            /*
             * W = A1 + V' * A2:
             *      W = A1
             *      W = W + V' * A2
             *
             */
            cudaMemcpy2DAsync( workW, sizeof(cuDoubleComplex) * ldW,
                               A1,    sizeof(cuDoubleComplex) * LDA1,
                               sizeof(cuDoubleComplex) * K, N1,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == ChamColumnwise ? ChamConjTrans : ChamNoTrans;
            transA2 = storev == ChamColumnwise ? ChamNoTrans : ChamConjTrans;

            rc = cublasZgemm( handle,
                              chameleon_cublas_const(transW), chameleon_cublas_const(ChamNoTrans),
                              K, N1, M2,
                              CUBLAS_SADDR(zone), workV /* M2*K  */, ldV,
                                                  A2    /* M2*N2 */, LDA2,
                              CUBLAS_SADDR(zone), workW /* K *N2 */, ldW );
            assert( rc == CUBLAS_STATUS_SUCCESS );
            info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

            if ( workC == NULL ) {
                /* W = op(T) * W */
                CUDA_ztrmm( ChamLeft, ChamUpper, trans, ChamNonUnit,
                            K, N2,
                            &zone, T,     LDT,
                                   workW, ldW,
                            handle );

                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    rc = cublasZaxpy( handle,
                                      K, CUBLAS_SADDR(mzone),
                                      workW + ldW  * j, 1,
                                      A1    + LDA1 * j, 1 );
                    assert( rc == CUBLAS_STATUS_SUCCESS );
                    info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
                }

                /* A2 = A2 - op(V) * W  */
                rc = cublasZgemm( handle,
                                  chameleon_cublas_const(transA2), chameleon_cublas_const(ChamNoTrans),
                                  M2, N2, K,
                                  CUBLAS_SADDR(mzone), workV /* M2 * K  */, ldV,
                                                       workW /* K  * N2 */, ldW,
                                  CUBLAS_SADDR(zone),  A2    /* M2 * N2 */, LDA2 );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

            } else {
                /* Wc = V * op(T) */
                rc = cublasZgemm( handle,
                                  chameleon_cublas_const(transA2), chameleon_cublas_const(trans),
                                  M2, K, K,
                                  CUBLAS_SADDR(zone),  workV, ldV,
                                                       T,     LDT,
                                  CUBLAS_SADDR(zzero), workC, ldC );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

                /* A1 = A1 - opt(T) * W */
                rc = cublasZgemm( handle,
                                  chameleon_cublas_const(trans), chameleon_cublas_const(ChamNoTrans),
                                  K, N1, K,
                                  CUBLAS_SADDR(mzone), T,     LDT,
                                                       workW, ldW,
                                  CUBLAS_SADDR(zone),  A1,    LDA1 );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

                /* A2 = A2 - Wc * W */
                rc = cublasZgemm( handle,
                                  chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(ChamNoTrans),
                                  M2, N2, K,
                                  CUBLAS_SADDR(mzone), workC, ldC,
                                                       workW, ldW,
                                  CUBLAS_SADDR(zone),  A2,    LDA2 );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
            }
        }
        else {
            /*
             * Column or Rowwise / Forward / Right
             * -----------------------------------
             *
             * Form  H * A  or  H' * A  where A  = ( A1 A2 )
             *
             */

            /*
             * Store in WORK (M1 == M2):
             *    - Workspace W for the copy of A1 + A2 * V' (M1 x K )
             *    - Workspace C for the copy of V * T        (K  x N2)
             *    - Workspace V for the copy of V            (K  x N2)
             */
            workW = WORK;
            ldW = M1;

            workC = workW + M1 * K;
            ldC = K;

            if ( L == 0 ) {
                workV = (cuDoubleComplex*)V;
                ldV   = LDV;
            }
            else {
                if ( LWORK < wrsize ) {
                    workC = NULL;
                    workV = workW + M2 * K;
                }
                else {
                    workV = workC + K * N2;
                }

                if ( storev == ChamColumnwise ) {
                    ldV = N2;

                    /*
                     * Backup V, and put 0 in the lower part
                     */
                    cudaMemcpy2DAsync( workV, sizeof(cuDoubleComplex) * ldV,
                                       V,     sizeof(cuDoubleComplex) * LDV,
                                       sizeof(cuDoubleComplex) * N2, K,
                                       cudaMemcpyDeviceToDevice, stream );

                    j = 0;
                    i = N2 - L + 1;
                    for(; (i < N2) && (j < K); i++, j++ ) {
                        cudaMemsetAsync( workV + j * ldV + i, 0,
                                         sizeof(cuDoubleComplex) * (N2 - i),
                                         stream );
                    }
                }
                else {
                    ldV = K;

                    /*
                     * Backup V, and put 0 in the upper part
                     */
                    cudaMemcpy2DAsync( workV, sizeof(cuDoubleComplex) * ldV,
                                       V,     sizeof(cuDoubleComplex) * LDV,
                                       sizeof(cuDoubleComplex) * K, N2,
                                       cudaMemcpyDeviceToDevice, stream );

                    for(j = 1; j < K; j++) {
                        cudaMemsetAsync( workV + ldV * ( N2 - L + j ), 0,
                                         sizeof(cuDoubleComplex) * j,
                                         stream );
                    }
                }
            }

            /*
             * W = A1 + A2 * V':
             *      W = A1
             *      W = W + A2 * V'
             *
             */
            cudaMemcpy2DAsync( workW, sizeof(cuDoubleComplex) * ldW,
                               A1,    sizeof(cuDoubleComplex) * LDA1,
                               sizeof(cuDoubleComplex) * M1, K,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == ChamColumnwise ? ChamNoTrans : ChamConjTrans;
            transA2 = storev == ChamColumnwise ? ChamConjTrans : ChamNoTrans;

            rc = cublasZgemm( handle,
                              chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(transW),
                              M1, K, N2,
                              CUBLAS_SADDR(zone), A2    /* M1*N2 */, LDA2,
                                                  workV /* K *N2 */, ldV,
                              CUBLAS_SADDR(zone), workW /* M1*K  */, ldW);
            assert( rc == CUBLAS_STATUS_SUCCESS );
            info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

            if ( workC == NULL ) {
                /* W = W * op(T) */
                CUDA_ztrmm( ChamRight, ChamUpper, trans, ChamNonUnit,
                            M2, K,
                            &zone, T,     LDT,
                                   workW, ldW,
                            handle );

                /* A1 = A1 - W = A1 - W * op(T) */
                for(j = 0; j < K; j++) {
                    rc = cublasZaxpy( handle,
                                      M1, CUBLAS_SADDR(mzone),
                                      workW + ldW  * j, 1,
                                      A1    + LDA1 * j, 1 );
                    assert( rc == CUBLAS_STATUS_SUCCESS );
                    info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
                }

                /* A2 = A2 - W * op(V)  */
                rc = cublasZgemm( handle,
                                  chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(transA2),
                                  M2, N2, K,
                                  CUBLAS_SADDR(mzone), workW /* M2*K  */, ldW,
                                                       workV /* K *N2 */, ldV,
                                  CUBLAS_SADDR(zone),  A2    /* M2*N2 */, LDA2);
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

            } else {
                /* A1 = A1 - W * opt(T) */
                rc = cublasZgemm( handle,
                                   chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(trans),
                                   M1, K, K,
                                   CUBLAS_SADDR(mzone), workW, ldW,
                                                        T,    LDT,
                                   CUBLAS_SADDR(zone),  A1,   LDA1 );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

                /* Wc = op(T) * V */
                rc = cublasZgemm( handle,
                                   chameleon_cublas_const(trans), chameleon_cublas_const(transA2),
                                   K, N2, K,
                                   CUBLAS_SADDR(zone),  T,     LDT,
                                                        workV, ldV,
                                   CUBLAS_SADDR(zzero), workC, ldC );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;

                /* A2 = A2 - W * Wc */
                rc = cublasZgemm( handle,
                                   chameleon_cublas_const(ChamNoTrans), chameleon_cublas_const(ChamNoTrans),
                                   M2, N2, K,
                                   CUBLAS_SADDR(mzone), workW, ldW,
                                                        workC, ldC,
                                   CUBLAS_SADDR(zone),  A2,    LDA2 );
                assert( rc == CUBLAS_STATUS_SUCCESS );
                info += (rc == CUBLAS_STATUS_SUCCESS) ? 0 : 1;
            }
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }

    return (info == 0) ? CHAMELEON_SUCCESS : CHAMELEON_ERR_UNEXPECTED;
}
