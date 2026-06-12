/**
 *
 * @file parsec/codelet_ztplqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztplqt PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztplqt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int M;
    int N;
    int L;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *WORK;

    parsec_dtd_unpack_args(
        this_task, &M, &N, &L, &ib, &A, &lda, &B, &ldb, &T, &ldt, &WORK );

    CORE_zlaset( ChamUpperLower, ib, M, 0., 0., T, ldt );
    CORE_ztplqt( M, N, L, ib,
                 A, lda, B, ldb, T, ldt, WORK );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztplqt_parsec, options->priority, "tplqt",
        sizeof(int),   &M,   VALUE,
        sizeof(int),   &N,   VALUE,
        sizeof(int),   &L,   VALUE,
        sizeof(int),   &ib,  VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,  RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PASSED_BY_REF,  RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | OUTPUT,
        sizeof(int), &(tileT->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*(ib+1)*nb, NULL, SCRATCH,
        PARSEC_DTD_ARG_END );
}
