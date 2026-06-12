/**
 *
 * @file quark/codelet_zsytrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Marc Sergent
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zsytrf_nopiv_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int n;
    CHAM_tile_t *tileA;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int iinfo;
    int info = 0;

    quark_unpack_args_6(quark, uplo, n, tileA, sequence, request, iinfo);
    info = TCORE_zsytf2_nopiv(uplo, n, tileA);
    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, iinfo+info );
    }
}

void INSERT_TASK_zsytrf_nopiv(const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_POTRF;
    QUARK_Insert_Task(opt->quark, CORE_zsytrf_nopiv_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,                VALUE,
        sizeof(int),                     &n,                   VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INOUT,
        sizeof(RUNTIME_sequence_t*),       &(options->sequence), VALUE,
        sizeof(RUNTIME_request_t*),        &(options->request),  VALUE,
        sizeof(int),                     &iinfo,               VALUE,
        0);

    (void)nb;
}
