/**
 *
 * @file pzplrnt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  chameleon_pzplrnt - Generate a random matrix by tiles.
 */
void chameleon_pzplrnt_generic( CHAM_desc_t            *A,
                                int                     bigM,
                                int                     m0,
                                int                     n0,
                                unsigned long long int  seed,
                                RUNTIME_sequence_t     *sequence,
                                RUNTIME_request_t      *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n;
    int tempmm, tempnn;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        for (n = 0; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            INSERT_TASK_zplrnt(
                &options,
                tempmm, tempnn, A(m, n),
                bigM, m*A->mb + m0, n*A->nb + n0, seed );
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}

void chameleon_pzplrnt( CHAM_desc_t            *A,
                        int                     bigM,
                        int                     m0,
                        int                     n0,
                        unsigned long long int  seed,
                        RUNTIME_sequence_t     *sequence,
                        RUNTIME_request_t      *request )
{
    chameleon_pzplrnt_generic( A, bigM, m0, n0, seed, sequence, request );

    /* Mark written data for synchronization */
    A->sync = 1;
}
