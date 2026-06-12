/**
 *
 * @file pzsytrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf parallel algorithm
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @author Samuel Thibault
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void chameleon_pzsytrf(cham_uplo_t uplo, CHAM_desc_t *A,
                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempmm, tempnn;
    size_t ws_host   = 0;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    RUNTIME_options_ws_alloc( &options, 0, ws_host );

    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );

            INSERT_TASK_zsytrf_nopiv(
                &options,
                ChamLower, tempkm, A->mb,
                A(k, k), A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, ChamLower, ChamTrans, ChamNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k),
                          A(m, k));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_zsyrk(
                    &options,
                    ChamLower, ChamNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k),
                     1.0, A(n, n));

                for (m = n+1; m < A->mt; m++) {
                    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k),
                               A(n, k),
                        zone,  A(m, n));
                }
                chameleon_data_flush( sequence, A(n, k), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /*
     *  ChamUpper
     */
    else {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = A->get_blkdim( A, k, DIM_n, A->n );
            INSERT_TASK_zsytrf_nopiv(
                &options,
                ChamUpper,
                tempkm, A->mb,
                A(k, k), A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, ChamUpper, ChamTrans, ChamNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k),
                          A(k, n));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                INSERT_TASK_zsyrk(
                    &options,
                    ChamUpper, ChamTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m),
                     1.0, A(m, m));

                for (n = m+1; n < A->nt; n++) {
                    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                    INSERT_TASK_zgemm(
                        &options,
                        ChamTrans, ChamNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m),
                               A(k, n),
                        zone,  A(m, n));
                }
                chameleon_data_flush( sequence, A(k, m), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
