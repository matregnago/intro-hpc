/**
 *
 * @file pzplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Mathis Rade
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n

/**
 * @brief Generic tile algorithm to generate an hermitian random matrix
 *
 * This is the version to use by default.
 */
void chameleon_pzplghe_generic( double                 bump,
                                cham_uplo_t            uplo,
                                CHAM_desc_t           *A,
                                int                    bigM,
                                int                    m0,
                                int                    n0,
                                unsigned long long int seed,
                                RUNTIME_sequence_t    *sequence,
                                RUNTIME_request_t     *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, minmn;
    int tempmm, tempnn;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    minmn = chameleon_min( A->mt, A->nt );
    switch ( uplo ) {
    case ChamLower:
        for (n = 0; n < minmn; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            for (m = n; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    bigM, m*A->mb + m0, n*A->nb + n0, seed );
            }
        }
        break;

    case ChamUpper:
        for (m = 0; m < minmn; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            for (n = m; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    bigM, m*A->mb + m0, n*A->nb + n0, seed );
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            for (n = 0; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    bigM, m*A->mb + m0, n*A->nb + n0, seed );
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}

void chameleon_pzplghe( double bump, cham_uplo_t uplo, CHAM_desc_t *A,
                        int bigM, int m0, int n0, unsigned long long int seed,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    chameleon_pzplghe_generic( bump, uplo, A, bigM, m0, n0, seed, sequence, request );
    /* Mark written data for synchronization */
    A->sync = 1;
}
