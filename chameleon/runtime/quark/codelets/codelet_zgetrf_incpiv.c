/**
 *
 * @file quark/codelet_zgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Raphael Boucherie
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgetrf_incpiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    int *IPIV;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    cham_bool_t check_info;
    int iinfo;

    int info;

    quark_unpack_args_9(quark, m, n, ib, tileA, IPIV, sequence, request, check_info, iinfo);
    TCORE_zgetrf_incpiv(m, n, ib, tileA, IPIV, &info);
    if ( (info != CHAMELEON_SUCCESS) && check_info ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, iinfo+info );
    }
}

void INSERT_TASK_zgetrf_incpiv(const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *L, int Lm, int Ln,
                              int *IPIV,
                              cham_bool_t check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(opt->quark, CORE_zgetrf_incpiv_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                     INOUT,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(RUNTIME_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(RUNTIME_request_t*),            &(options->request),       VALUE,
        sizeof(cham_bool_t),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);

    (void)L;
    (void)Lm;
    (void)Ln;
    (void)nb;
}
