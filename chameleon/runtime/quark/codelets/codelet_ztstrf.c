/**
 *
 * @file quark/codelet_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf Quark codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"
#include "coreblas/cblas.h"
#include <math.h>

void CORE_ztstrf_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    int nb;
    CHAM_tile_t *tileU;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileL;
    int *IPIV;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    cham_bool_t check_info;
    int iinfo;

    int info;

    quark_unpack_args_14(quark, m, n, ib, nb, tileU, tileA, tileL, IPIV, WORK, ldwork, sequence, request, check_info, iinfo);
    TCORE_ztstrf(m, n, ib, nb, tileU, tileA, tileL, IPIV, WORK, ldwork, &info);
    if ( (info != CHAMELEON_SUCCESS) && check_info ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, iinfo+info );
    }
}

void INSERT_TASK_ztstrf(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *U, int Um, int Un,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *L, int Lm, int Ln,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TSTRF;
    QUARK_Insert_Task(opt->quark, CORE_ztstrf_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(int),                        &nb,            VALUE,
        sizeof(void*), RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),                     INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                     INOUT | LOCALITY,
        sizeof(void*), RTBLKADDR(L, CHAMELEON_Complex64_t, Lm, Ln),                     OUTPUT,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,    NULL,                  SCRATCH,
        sizeof(int),                        &nb,            VALUE,
        sizeof(RUNTIME_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(RUNTIME_request_t*),            &(options->request),       VALUE,
        sizeof(cham_bool_t),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
