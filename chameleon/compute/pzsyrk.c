/**
 *
 * @file pzsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
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
#define C(m,n) C,  m,  n
/**
 *  Parallel tile symmetric rank-k update - dynamic scheduling
 */
void chameleon_pzsyrk(cham_uplo_t uplo, cham_trans_t trans,
                         CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                         CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, k;
    int tempnn, tempmm, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for (n = 0; n < C->nt; n++) {
        tempnn = C->get_blkdim( C, n, DIM_n, C->n );
        /*
         *  ChamNoTrans
         */
        if (trans == ChamNoTrans) {
            for (k = 0; k < A->nt; k++) {
                tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                zbeta = k == 0 ? beta : zone;
                INSERT_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkn, A->mb,
                    alpha, A(n, k), /* ldan * K */
                    zbeta, C(n, n)); /* ldc  * N */
            }
            /*
             *  ChamNoTrans / ChamLower
             */
            if (uplo == ChamLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                    for (k = 0; k < A->nt; k++) {
                        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* ldam * K */
                                   A(n, k),  /* ldan * K */
                            zbeta, C(m, n)); /* ldc  * N */
                    }
                }
            }
            /*
             *  ChamNoTrans / ChamUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                    for (k = 0; k < A->nt; k++) {
                        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamTrans,
                            tempnn, tempmm, tempkn, A->mb,
                            alpha, A(n, k),  /* ldan * K */
                                   A(m, k),  /* ldam * M */
                            zbeta, C(n, m)); /* ldc  * M */
                    }
                }
            }
        }
        /*
         *  ChamTrans
         */
        else {
            for (k = 0; k < A->mt; k++) {
                tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                zbeta = k == 0 ? beta : zone;
                INSERT_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkm, A->mb,
                    alpha, A(k, n),  /* lda * N */
                    zbeta, C(n, n)); /* ldc * N */
            }
            /*
             *  ChamTrans / ChamLower
             */
            if (uplo == ChamLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                    for (k = 0; k < A->mt; k++) {
                        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamNoTrans,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * M */
                                   A(k, n),  /* lda * N */
                            zbeta, C(m, n)); /* ldc * N */
                    }
                }
            }
            /*
             *  ChamTrans / ChamUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );
                    for (k = 0; k < A->mt; k++) {
                        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamNoTrans,
                            tempnn, tempmm, tempkm, A->mb,
                            alpha, A(k, n),  /* lda * K */
                                   A(k, m),  /* lda * M */
                            zbeta, C(n, m)); /* ldc * M */
                    }
                }
            }
        }
    }
    /* Mark written data for synchronization */
    C->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
