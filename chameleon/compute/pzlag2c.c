/**
 *
 * @file pzlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions mixed zc -> ds
 *
 */
#include "control/common.h"

#define A(  _m_, _n_ ) A,  (_m_), (_n_)
#define B(  _m_, _n_ ) B,  (_m_), (_n_)

/**
 *
 */
void chameleon_pclag2z( CHAM_desc_t *A, CHAM_desc_t *B,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int tempmm, tempnn;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for(m = 0; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        for(n = 0; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            INSERT_TASK_clag2z(
                &options,
                tempmm, tempnn, A->mb,
                A(m, n),
                B(m, n));
        }
    }
    /* Mark written data for synchronization */
    B->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}

/**
 *
 */
void chameleon_pzlag2c( CHAM_desc_t *A, CHAM_desc_t *B,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int tempmm, tempnn;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for(m = 0; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        for(n = 0; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            INSERT_TASK_zlag2c(
                &options,
                tempmm, tempnn, A->mb,
                A(m, n),
                B(m, n));
        }
    }
    /* Mark written data for synchronization */
    B->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
