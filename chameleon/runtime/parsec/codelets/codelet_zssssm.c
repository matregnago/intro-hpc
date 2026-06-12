/**
 *
 * @file parsec/codelet_zssssm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zssssm PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zssssm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
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
    CHAMELEON_Complex64_t *L1;
    int ldl1;
    CHAMELEON_Complex64_t *L2;
    int ldl2;
    int *IPIV;

    parsec_dtd_unpack_args(
        this_task, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &L1, &ldl1, &L2, &ldl2, &IPIV );

    CORE_zssssm( m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zssssm(const RUNTIME_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n,
                       const CHAM_desc_t *A2, int A2m, int A2n,
                       const CHAM_desc_t *L1, int L1m, int L1n,
                       const CHAM_desc_t *L2, int L2m, int L2n,
                       const int *IPIV)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA1 = A1->get_blktile( A1, A1m, A1n );
    CHAM_tile_t *tileA2 = A2->get_blktile( A2, A2m, A2n );
    CHAM_tile_t *tileL1 = L1->get_blktile( L1, L1m, L1n );
    CHAM_tile_t *tileL2 = L2->get_blktile( L2, L2m, L2n );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zssssm_parsec, options->priority, "ssssm",
        sizeof(int),           &m1,                                VALUE,
        sizeof(int),           &n1,                                VALUE,
        sizeof(int),           &m2,                                VALUE,
        sizeof(int),           &n2,                                VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A1, CHAMELEON_Complex64_t, A1m, A1n ), chameleon_parsec_get_arena_index( A1 ) | INOUT,
        sizeof(int), &(tileA1->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( A2, CHAMELEON_Complex64_t, A2m, A2n ), chameleon_parsec_get_arena_index( A2 ) | INOUT | AFFINITY,
        sizeof(int), &(tileA2->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( L1, CHAMELEON_Complex64_t, L1m, L1n ), chameleon_parsec_get_arena_index( L1 ) | INPUT,
        sizeof(int), &(tileL1->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( L2, CHAMELEON_Complex64_t, L2m, L2n ), chameleon_parsec_get_arena_index( L2 ) | INPUT,
        sizeof(int), &(tileL2->ld), VALUE,
        sizeof(int*),          &IPIV,                              VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
