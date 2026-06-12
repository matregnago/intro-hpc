/**
 *
 * @file pzlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n

void chameleon_pzlacpy_generic( cham_uplo_t         uplo,
                                CHAM_desc_t        *A,
                                CHAM_desc_t        *B,
                                RUNTIME_sequence_t *sequence,
                                RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int X, Y;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch (uplo) {
    /*
     *  ChamUpper
     */
    case ChamUpper:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            if (m < A->nt) {
                Y = m == A->nt-1 ? A->n-m*A->nb : A->nb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpper,
                    X, Y,
                    A(m, m),
                    B(m, m));
            }
            for (n = m+1; n < A->nt; n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower,
                    X, Y,
                    A(m, n),
                    B(m, n));
            }
        }
        break;
    /*
     *  ChamLower
     */
    case ChamLower:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            if (m < A->nt) {
                Y = m == A->nt-1 ? A->n-m*A->nb : A->nb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamLower,
                    X, Y,
                    A(m, m),
                    B(m, m));
            }
            for (n = 0; n < chameleon_min(m, A->nt); n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower,
                    X, Y,
                    A(m, n),
                    B(m, n));
            }
        }
        break;
    /*
     *  ChamUpperLower
     */
    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            for (n = 0; n < A->nt; n++) {
                Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower,
                    X, Y,
                    A(m, n),
                    B(m, n));
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}

void chameleon_pzlacpy( cham_uplo_t         uplo,
                        CHAM_desc_t        *A,
                        CHAM_desc_t        *B,
                        RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t  *request )
{
    chameleon_pzlacpy_generic( uplo, A, B, sequence, request );
    /* Mark written data for synchronization */
    B->sync = 1;
}
