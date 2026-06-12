/**
 *
 * @file pztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 */
void chameleon_pztrtri(cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkn, tempkm, tempmm, tempnn;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);
    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, diag,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k),
                           A(m, k));
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                for (n = 0; n < k; n++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k),
                              A(k, n),
                        zone, A(m, n));
                }
                chameleon_data_flush( sequence, A(m, k), request->flush );
            }
            for (n = 0; n < k; n++) {
                chameleon_data_flush( sequence, A(k, n), request->flush );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, diag,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k),
                          A(k, n));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_ztrtri(
                &options,
                uplo, diag,
                tempkn, A->mb,
                A(k, k), A->nb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /*
     *  ChamUpper
     */
    else {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, diag,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k),
                           A(k, n));
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                for (m = 0; m < k; m++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k),
                              A(k, n),
                        zone, A(m, n));
                }
                chameleon_data_flush( sequence, A(k, n), request->flush );
            }
            for (m = 0; m < k; m++) {
                chameleon_data_flush( sequence, A(m, k), request->flush );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, diag,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k),
                          A(m, k));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_ztrtri(
                &options,
                uplo, diag,
                tempkm, A->mb,
                A(k, k), A->mb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
