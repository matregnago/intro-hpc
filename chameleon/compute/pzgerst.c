/**
 *
 * @file pzgerst.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgerst parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> d
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n

void chameleon_pzgerst( cham_uplo_t         uplo,
                        CHAM_desc_t        *A,
                        RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for(m = 0; m < A->mt; m++) {
        int tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        int nmin   = ( uplo == ChamUpper ) ? m                         : 0;
        int nmax   = ( uplo == ChamLower ) ? chameleon_min(m+1, A->nt) : A->nt;

        for(n = nmin; n < nmax; n++) {
            CHAM_tile_t *tile = A->get_blktile( A, m, n );

            if (( tile->rank == A->myrank ) &&
                ( tile->flttype != ChamComplexDouble ) )
            {
                int tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                INSERT_TASK_zgerst( &options,
                                     tempmm, tempnn, A( m, n ) );
            }
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
