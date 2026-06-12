/**
 *
 * @file pzlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal parallel algorithm
 *
 * @version 1.4.0
 * @author Dalal Sukkari
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  Parallel scale of a matrix A
 */
void chameleon_pzlascal(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int minmnt = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch(uplo) {
    case ChamLower:
        for (n = 0; n < minmnt; n++) {
            tempnm = A->get_blkdim( A, n, DIM_m, A->m );
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            INSERT_TASK_zlascal(
                &options,
                ChamLower, tempnm, tempnn, A->mb,
                alpha, A(n, n));

            for (m = n+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n));
            }
        }
        break;

    case ChamUpper:
        for (m = 0; m < minmnt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            tempmn = A->get_blkdim( A, m, DIM_n, A->n );

            INSERT_TASK_zlascal(
                &options,
                ChamUpper, tempmm, tempmn, A->mb,
                alpha, A(m, m));

            for (n = m+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n));
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            for (n = 0; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                INSERT_TASK_zlascal(
                    &options,
                    ChamUpperLower, tempmm, tempnn, A->mb,
                    alpha, A(m, n));
            }
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
