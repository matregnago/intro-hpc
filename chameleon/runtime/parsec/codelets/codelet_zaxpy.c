/**
 *
 * @file parsec/codelet_zaxpy.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy PaRSEC codelet
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
CORE_zaxpy_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    int M;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int incA;
    CHAMELEON_Complex64_t *B;
    int incB;

    parsec_dtd_unpack_args(
        this_task, &M, &alpha, &A, &incA, &B, &incB );

    CORE_zaxpy( M, alpha, A, incA, B, incB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zaxpy(const RUNTIME_option_t *options,
                      int M, CHAMELEON_Complex64_t alpha,
                      const CHAM_desc_t *A, int Am, int An, int incA,
                      const CHAM_desc_t *B, int Bm, int Bn, int incB)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zaxpy_parsec, options->priority,  "axpy",
        sizeof(int),                   &M,     VALUE,
        sizeof(CHAMELEON_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF, RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),                   &incA,  VALUE,
        PASSED_BY_REF, RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int),                   &incB,  VALUE,
        PARSEC_DTD_ARG_END );
}
