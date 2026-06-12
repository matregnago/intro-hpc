/**
 *
 * @file pzher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile Hermitian rank-k update - dynamic scheduling
 */
void chameleon_pzher2k( cham_uplo_t uplo, cham_trans_t trans,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        double beta,  CHAM_desc_t *C,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, k, mmin, mmax;
    int tempnn, tempmm, tempkn, tempkm;

    CHAMELEON_Complex64_t zone   = (CHAMELEON_Complex64_t)1.0;
    CHAMELEON_Complex64_t zbeta;
    double dbeta;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for (n = 0; n < C->nt; n++) {
        tempnn = C->get_blkdim( C, n, DIM_n, C->n );

        if (uplo == ChamLower) {
            mmin = n+1;
            mmax = C->mt;
        }
        else {
            mmin = 0;
            mmax = n;
        }

        /*
         *  ChamNoTrans
         */
        if (trans == ChamNoTrans) {
            for (k = 0; k < A->nt; k++) {
                tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                dbeta = k == 0 ? beta : 1.0;
                INSERT_TASK_zher2k(
                    &options,
                    uplo, trans,
                    tempnn, tempkn, A->mb,
                    alpha, A(n, k), /* ldan * K */
                           B(n, k),
                    dbeta, C(n, n)); /* ldc  * N */
            }
            for (m = mmin; m < mmax; m++) {
                tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                for (k = 0; k < A->nt; k++) {
                    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                    zbeta = k == 0 ? (CHAMELEON_Complex64_t)beta : zone;
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        tempmm, tempnn, tempkn, A->mb,
                        alpha, A(m, k),
                               B(n, k),
                        zbeta, C(m, n));

                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        tempmm, tempnn, tempkn, A->mb,
                        conj(alpha), B(m, k),
                                     A(n, k),
                        zone,        C(m, n));
                }
            }
        }
        /*
         *  ChamConjTrans
         */
        else {
            for (k = 0; k < A->mt; k++) {
                tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                dbeta = k == 0 ? beta : 1.0;
                INSERT_TASK_zher2k(
                    &options,
                    uplo, trans,
                    tempnn, tempkm, A->mb,
                    alpha, A(k, n),  /* lda * N */
                           B(k, n),
                    dbeta, C(n, n)); /* ldc * N */
            }
            for (m = mmin; m < mmax; m++) {
                tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                for (k = 0; k < A->mt; k++) {
                    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                    zbeta = k == 0 ? (CHAMELEON_Complex64_t)beta : zone;
                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        tempmm, tempnn, tempkm, A->mb,
                        alpha, A(k, m),
                               B(k, n),
                        zbeta, C(m, n));

                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        tempmm, tempnn, tempkm, A->mb,
                        conj(alpha), B(k, m),
                                     A(k, n),
                        zone,        C(m, n) );
                }
            }
        }
    }
    /* Mark written data for synchronization */
    C->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
