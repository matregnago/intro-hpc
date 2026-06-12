/**
 *
 * @file parsec/codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_blas PaRSEC wrapper
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zherfb_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int n;
    int k;
    int ib;
    int nb;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task,   &uplo,   &n,   &k,   &ib,   &nb, &A,   &lda, &T,   &ldt, &C,   &ldc, &WORK,   &ldwork);

    CORE_zherfb( uplo, n, k, ib, nb,
                A, lda, T, ldt,
                C, ldc, WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zherfb_parsec, options->priority, "herfb",
        sizeof(cham_uplo_t), &uplo, VALUE,
        sizeof(int),        &n,    VALUE,
        sizeof(int),        &k,    VALUE,
        sizeof(int),        &ib,   VALUE,
        sizeof(int),        &nb,   VALUE,
        PASSED_BY_REF,       RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INOUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn), INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*2*nb*nb,  NULL, SCRATCH,
        sizeof(int),        &nb,   VALUE,
        PARSEC_DTD_ARG_END );
}
