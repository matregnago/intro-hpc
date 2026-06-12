/**
 *
 * @file pzlatro.c
 *
 * @copyright 2026-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Matrix transposition parallel algorithm
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2026-01-23
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A, m, n
#define B(m,n) B, m, n

static inline void
chameleon_pzlatro_upperlower( CHAM_context_t *chamctxt,
                              cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *B,
                              RUNTIME_option_t *options )
{
    int m, n, tempm, tempn;

    for (n = 0; n < A->nt; n++) {
        tempn = A->get_blkdim( A, n, DIM_n, A->n );
        for (m = 0; m < A->mt; m++) {
            tempm = A->get_blkdim( A, m, DIM_m, A->m );

            INSERT_TASK_zlatro(
                options, ChamUpperLower, trans,
                tempm, tempn, A->mb,
                A( m, n ), B( n, m ) );
        }
    }
    (void)chamctxt;
}

static inline void
chameleon_pzlatro_upper( CHAM_context_t *chamctxt,
                         cham_trans_t trans,
                         CHAM_desc_t *A, CHAM_desc_t *B,
                         RUNTIME_option_t *options )
{
    int m, n, tempm, tempn;

    for( n = 0; n < A->nt; n++) {
        tempn = A->get_blkdim( A, n, DIM_n, A->n );
        for( m = 0; m < chameleon_min( n + 1, A->mt ); m++) {
            tempm = A->get_blkdim( A, m, DIM_m, A->m );

            cham_uplo_t uplo = m == n ? ChamUpper : ChamUpperLower;

            INSERT_TASK_zlatro(
                options, uplo, trans,
                tempm, tempn, A->mb,
                A( m, n ), B( n, m ) );
        }
    }
    (void)chamctxt;
}

static inline void
chameleon_pzlatro_lower( CHAM_context_t *chamctxt,
                         cham_trans_t trans,
                         CHAM_desc_t *A, CHAM_desc_t *B,
                         RUNTIME_option_t *options )
{
    int m, n, tempm, tempn;

    for( n = 0; n < A->nt; n++) {
        tempn = A->get_blkdim( A, n, DIM_n, A->n );
        for( m = n; m < A->mt; m++) {
            tempm = A->get_blkdim( A, m, DIM_m, A->m );

            cham_uplo_t uplo = m == n ? ChamLower : ChamUpperLower;

            INSERT_TASK_zlatro(
                options, uplo, trans,
                tempm, tempn, A->mb,
                A( m, n ), B( n, m ) );
        }
    }
    (void)chamctxt;
}

/**
 *  Parallel tile transposition - dynamic scheduling
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = ChamUpper: the upper triangle of A and the lower triangle of B
 *          are referenced.
 *          = ChamLower: the lower triangle of A and the upper triangle of B
 *          are referenced.
 *          = ChamUpperLower: All A and B are referenced.
 */
void chameleon_pzlatro( cham_uplo_t uplo, cham_trans_t trans,
                        CHAM_desc_t *A, CHAM_desc_t *B,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch( uplo ) {
        case ChamUpperLower:
            chameleon_pzlatro_upperlower( chamctxt, trans, A, B, &options );
            break;
        case ChamUpper:
            chameleon_pzlatro_upper( chamctxt, trans, A, B, &options );
            break;
        case ChamLower:
            chameleon_pzlatro_lower( chamctxt, trans, A, B, &options );
            break;
    }
    /* Mark written data for synchronization */
    B->sync = 1;

    RUNTIME_options_finalize( &options, chamctxt );
}
