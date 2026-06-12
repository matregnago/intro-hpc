/**
 *
 * @file pzlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
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
/**
 *  Parallel initialization a 2-D array A to BETA on the diagonal and
 *  ALPHA on the offdiagonals.
 */
void chameleon_pzlaset( cham_uplo_t uplo,
                        CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                        CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int i, j;
    int tempim;
    int tempjm, tempjn;
    int minmn = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    if (uplo == ChamLower) {
       for (j = 0; j < minmn; j++){
           tempjm = A->get_blkdim( A, j, DIM_m, A->m );
           tempjn = A->get_blkdim( A, j, DIM_n, A->n );
           INSERT_TASK_zlaset(
               &options,
               ChamLower, tempjm, tempjn, alpha, beta,
               A(j, j));

           for (i = j+1; i < A->mt; i++){
               tempim = A->get_blkdim( A, i, DIM_m, A->m );
               INSERT_TASK_zlaset(
                   &options,
                   ChamUpperLower, tempim, tempjn, alpha, alpha,
                   A(i, j));
           }
       }
    }
    else if (uplo == ChamUpper) {
        for (i = 0; i < A->mt; i++) {
            tempim = A->get_blkdim( A, i, DIM_m, A->m );

            if ( i < A->nt ) {
                j = i;
                tempjn = A->get_blkdim( A, j, DIM_n, A->n );

                INSERT_TASK_zlaset(
                    &options,
                    uplo, tempim, tempjn,
                    alpha, beta, A(i, j));
            }
            for (j = i+1; j < A->nt; j++) {
                tempjn = A->get_blkdim( A, j, DIM_n, A->n );

                INSERT_TASK_zlaset(
                    &options,
                    ChamUpperLower, tempim, tempjn,
                    alpha, alpha, A(i, j));
            }
        }
    }
    else {
       for (i = 0; i < A->mt; i++){
           tempim = A->get_blkdim( A, i, DIM_m, A->m );
           for (j = 0; j < A->nt; j++){
               tempjn = A->get_blkdim( A, j, DIM_n, A->n );
               INSERT_TASK_zlaset(
                   &options,
                   ChamUpperLower, tempim, tempjn,
                   alpha, (i == j) ? beta : alpha,
                   A(i, j));
           }
       }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
