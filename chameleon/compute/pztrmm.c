/**
 *
 * @file pztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n



/**
 *  Parallel tile triangular matrix-matrix multiplication - dynamic scheduling
 */
void chameleon_pztrmm(cham_side_t side, cham_uplo_t uplo,
                         cham_trans_t trans, cham_diag_t diag,
                         CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;

    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);
    /*
     *  ChamLeft / ChamUpper / ChamNoTrans
     */
    if (side == ChamLeft) {
        if (uplo == ChamUpper) {
            if (trans == ChamNoTrans) {
                for (m = 0; m < B->mt; m++) {
                    tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                    for (n = 0; n < B->nt; n++) {
                        tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = m+1; k < A->mt; k++) {
                            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                            INSERT_TASK_zgemm(
                                &options,
                                trans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, A(m, k),
                                       B(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
            /*
             *  ChamLeft / ChamUpper / Cham[Conj]Trans
             */
            else {
                for (m = B->mt-1; m > -1; m--) {
                    tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                    for (n = 0; n < B->nt; n++) {
                        tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            INSERT_TASK_zgemm(
                                &options,
                                trans, ChamNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, A(k, m),
                                       B(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
        }
        /*
         *  ChamLeft / ChamLower / ChamNoTrans
         */
        else {
            if (trans == ChamNoTrans) {
                for (m = B->mt-1; m > -1; m--) {
                    tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                    for (n = 0; n < B->nt; n++) {
                        tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = 0; k < m; k++) {
                            INSERT_TASK_zgemm(
                                &options,
                                trans, ChamNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, A(m, k),
                                       B(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
            /*
             *  ChamLeft / ChamLower / Cham[Conj]Trans
             */
            else {
                for (m = 0; m < B->mt; m++) {
                    tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                    for (n = 0; n < B->nt; n++) {
                        tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(m, m),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = m+1; k < A->mt; k++) {
                            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                            INSERT_TASK_zgemm(
                                &options,
                                trans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m),
                                       B(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
        }
    }
    /*
     *  ChamRight / ChamUpper / ChamNoTrans
     */
    else {
        if (uplo == ChamUpper) {
            if (trans == ChamNoTrans) {
                for (n = B->nt-1; n > -1; n--) {
                    tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                    for (m = 0; m < B->mt; m++) {
                        tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, B(m, k),
                                       A(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
            /*
             *  ChamRight / ChamUpper / Cham[Conj]Trans
             */
            else {
                for (n = 0; n < B->nt; n++) {
                    tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                    for (m = 0; m < B->mt; m++) {
                        tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = n+1; k < A->mt; k++) {
                            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, trans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),
                                       A(n, k),
                                zone,  B(m, n));
                        }
                    }
                }
            }
        }
        /*
         *  ChamRight / ChamLower / ChamNoTrans
         */
        else {
            if (trans == ChamNoTrans) {
                for (n = 0; n < B->nt; n++) {
                    tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                    for (m = 0; m < B->mt; m++) {
                        tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = n+1; k < A->mt; k++) {
                            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, trans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),
                                       A(k, n),
                                zone,  B(m, n));
                        }
                    }
                }
            }
            /*
             *  ChamRight / ChamLower / Cham[Conj]Trans
             */
            else {
                for (n = B->nt-1; n > -1; n--) {
                    tempnn = B->get_blkdim( B, n, DIM_n, B->n );
                    for (m = 0; m < B->mt; m++) {
                        tempmm = B->get_blkdim( B, m, DIM_m, B->m );
                        INSERT_TASK_ztrmm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempnn, A->mb,
                            alpha, A(n, n),  /* lda * tempkm */
                                   B(m, n)); /* ldb * tempnn */

                        for (k = 0; k < n; k++) {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                alpha, B(m, k),
                                       A(n, k),
                                zone,  B(m, n));
                        }
                    }
                }
            }
        }
    }
    /* Mark written data for synchronization */
    B->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
