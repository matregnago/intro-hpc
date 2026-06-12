/**
 *
 * @file pztile2band.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztile2band parallel algorithm
 *
 * @version 1.4.0
 * @author Azzam Haidar
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A, m, n
#define B(m,n)  B, m, n

/**
 *  Parallel copy of a band matrix from full NxN tile storage to band storage (LDABxN).
 */
void chameleon_pztile2band( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_tile_t *tileA;

    int k;
    int tempkm, tempkn;
    int minmnt = chameleon_min(A->mt, A->nt);
    int Bmb = B->mb;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /*
     *  ChamLower => Lower Band
     */
    if ( uplo == ChamLower ) {
        for (k = 0; k < minmnt; k++){
            tileA = A->get_blktile( A, k, k );

            /* Computes dimension on N with B since it is dimensioned with chameleon_min(A->m, A->n) */
            assert( A->i == B->j );
            assert( A->j >= B->j );

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            tempkn = B->get_blkdim( B, k, DIM_n, B->n );

            INSERT_TASK_zlaset( &options, ChamUpperLower, Bmb, tempkn,
                                0., 0., B, 0, k );

            INSERT_TASK_zlacpyx( &options, ChamLower, tempkm, tempkn,
                                 0, A, k, k, tileA->ld,
                                 0, B, 0, k, Bmb-1 );

            if ( k < minmnt-1 ) {
                tileA = A->get_blktile( A, k+1, k );

                tempkm = A->get_blkdim( A, k+1, DIM_m, A->m );

                INSERT_TASK_zlacpyx( &options, ChamUpper, tempkm, tempkn,
                                     0,     A, k+1, k, tileA->ld,
                                     Bmb-1, B, 0,   k, Bmb-1 );
            }
        }
    }
    else if ( uplo == ChamUpper ) {
        for (k = 0; k < minmnt; k++){
            tileA = A->get_blktile( A, k, k );

            /* Compute dimension on M with B since it is dimensioned with chameleon_min(A->m, A->n) */
            assert( A->i == B->i );
            assert( A->i >= B->j );

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            tempkn = B->get_blkdim( B, k, DIM_n, B->n );

            INSERT_TASK_zlaset( &options, ChamUpperLower, Bmb, tempkn,
                                0., 0., B, 0, k );

            INSERT_TASK_zlacpyx( &options, ChamUpper, tempkm, tempkn,
                                 0,     A, k, k, tileA->ld,
                                 Bmb-1, B, 0, k, Bmb-1 );

            if ( k > 0 ) {
                tileA = A->get_blktile( A, k-1, k );

                tempkm = A->get_blkdim( A, k-1, DIM_m, A->m );

                INSERT_TASK_zlacpyx( &options, ChamLower, tempkm, tempkn,
                                     0, A, k-1, k, tileA->ld,
                                     0, B, 0,   k, Bmb-1 );
            }
        }
    }
    /* Mark written data for synchronization */
    B->sync = 1;
    RUNTIME_options_finalize(&options, chamctxt);
}
#undef B
#undef A
