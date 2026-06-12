/**
 *
 * @file pmap.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-03-14
 *
 */
#include "control/common.h"

#define A(m, n) A, m, n

/**
 *  chameleon_pmap
 */
void chameleon_pmap( cham_uplo_t uplo, int ndata, cham_map_data_t *data,
                     cham_map_operator_t *op_fct, void *op_args,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    const CHAM_desc_t *A = data[0].desc;
    int i, m, n, minmn;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    minmn = chameleon_min( A->mt, A->nt );

    switch( uplo ) {
    case ChamUpper:
        for (m = 0; m < minmn; m++) {
            INSERT_TASK_map(
                &options, uplo, m, m,
                ndata, data,
                op_fct, op_args );

            for (n = m+1; n < A->nt; n++) {
                INSERT_TASK_map(
                    &options, ChamUpperLower, m, n,
                    ndata, data,
                    op_fct, op_args );
            }
        }
        break;

    case ChamLower:
        for (n = 0; n < minmn; n++){
            INSERT_TASK_map(
                &options, uplo, n, n,
                ndata, data,
                op_fct, op_args );
            for (m = n+1; m < A->mt; m++) {
                INSERT_TASK_map(
                    &options, ChamUpperLower, m, n,
                    ndata, data,
                    op_fct, op_args );
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            for (n = 0; n < A->nt; n++) {
                INSERT_TASK_map(
                    &options, ChamUpperLower, m, n,
                    ndata, data,
                    op_fct, op_args );
            }
        }
    }

    for( i=0; i<ndata; i++ ) {
        if ( data[i].access & ChamW ) {
            data[i].desc->sync = 1;
        }
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
