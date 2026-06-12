/**
 *
 * @file parsec/codelet_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmqr_hetra1 PaRSEC codelet
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztsmqr_hetra1_parsec( parsec_execution_stream_t *context,
                           parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A1;
    int lda1;
    CHAMELEON_Complex64_t *A2;
    int lda2;
    CHAMELEON_Complex64_t *V;
    int ldv;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &V, &ldv, &T, &ldt, &WORK, &ldwork);

    CORE_ztsmqr_hetra1( side, trans, m1, n1, m2, n2, k, ib,
                        A1, lda1, A2, lda2,
                        V, ldv, T, ldt,
                        WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztsmqr_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n,
                              const CHAM_desc_t *A2, int A2m, int A2n,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA1 = A1->get_blktile( A1, A1m, A1n );
    CHAM_tile_t *tileA2 = A2->get_blktile( A2, A2m, A2n );
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    int ldwork = side == ChamLeft ? ib : nb;

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztsmqr_hetra1_parsec, options->priority, "tsmqr_hetra1",
        sizeof(cham_side_t), &side,   VALUE,
        sizeof(cham_trans_t), &trans,  VALUE,
        sizeof(int),        &m1,     VALUE,
        sizeof(int),        &n1,     VALUE,
        sizeof(int),        &m2,     VALUE,
        sizeof(int),        &n2,     VALUE,
        sizeof(int),        &k,      VALUE,
        sizeof(int),        &ib,     VALUE,
        PASSED_BY_REF,       RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n), INOUT,
        sizeof(int), &(tileA1->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n), INOUT | AFFINITY,
        sizeof(int), &(tileA2->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(V,  CHAMELEON_Complex64_t, Vm,  Vn),  INPUT,
        sizeof(int), &(tileV->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(T,  CHAMELEON_Complex64_t, Tm,  Tn),  INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL, SCRATCH,
        sizeof(int),        &ldwork, VALUE,
        PARSEC_DTD_ARG_END );
}
