/**
 *
 * @file parsec/codelet_zunmqr.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr PaRSEC codelet
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

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "gpucublas.h"
#endif

static inline int
CORE_zunmqr_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m, &n, &k, &ib, &A, &lda, &T, &ldt, &C, &ldc, &WORK, &ldwork );

    CORE_zunmqr( side, trans, m, n, k, ib,
                 A, lda, T, ldt, C, ldc, WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
static __thread cublasHandle_t   cham_parsec_unmqr_cublas_handle = NULL;
static __thread cuDoubleComplex *cham_parsec_unmqr_scratch = NULL;
static __thread int              cham_parsec_unmqr_scratch_elts = 0;

static int
CORE_zunmqr_parsec_cuda( parsec_execution_stream_t *es,
                         parsec_task_t             *this_task,
                         void                      *cuda_stream )
{
    cham_side_t side;
    cham_trans_t trans;
    int m, n, k, ib, lda, ldt, ldc, ldwork;
    CHAMELEON_Complex64_t *A, *T, *C, *WORK_unused;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m, &n, &k, &ib,
        &A, &lda, &T, &ldt, &C, &ldc, &WORK_unused, &ldwork );

    A = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[0].data_out);
    T = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[1].data_out);
    C = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[2].data_out);

    if( NULL == cham_parsec_unmqr_cublas_handle ) {
        cublasStatus_t cs = cublasCreate( &cham_parsec_unmqr_cublas_handle );
        if( cs != CUBLAS_STATUS_SUCCESS ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
    }
    cublasSetStream( cham_parsec_unmqr_cublas_handle, (cudaStream_t)cuda_stream );

    /* CUDA_zunmqrt needs ib*ldwork workspace (StarPU sizes it ib*nb). */
    int need = ib * ldwork;
    if( need > cham_parsec_unmqr_scratch_elts ) {
        if( cham_parsec_unmqr_scratch != NULL ) {
            cudaFree( cham_parsec_unmqr_scratch );
        }
        if( cudaMalloc( (void **)&cham_parsec_unmqr_scratch,
                        (size_t)need * sizeof(cuDoubleComplex) ) != cudaSuccess ) {
            cham_parsec_unmqr_scratch = NULL;
            cham_parsec_unmqr_scratch_elts = 0;
            return PARSEC_HOOK_RETURN_ERROR;
        }
        cham_parsec_unmqr_scratch_elts = need;
    }

    CUDA_zunmqrt( side, trans, m, n, k, ib,
                  (const cuDoubleComplex *)A, lda,
                  (const cuDoubleComplex *)T, ldt,
                  (cuDoubleComplex *)C, ldc,
                  cham_parsec_unmqr_scratch, ldwork,
                  cham_parsec_unmqr_cublas_handle );

    (void)es;
    (void)WORK_unused;
    return PARSEC_HOOK_RETURN_DONE;
}
#endif /* CHAMELEON_USE_CUDA && !CHAMELEON_SIMULATION */

void INSERT_TASK_zunmqr(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m, int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    parsec_dtd_taskpool_insert_task_cuda(
        PARSEC_dtd_taskpool,
        CORE_zunmqr_parsec, CORE_zunmqr_parsec_cuda,
        options->priority, "unmqr",
        sizeof(cham_side_t),    &side,                              VALUE,
        sizeof(cham_trans_t),    &trans,                             VALUE,
        sizeof(int),           &m,                                 VALUE,
        sizeof(int),           &n,                                 VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,   NULL,                          SCRATCH,
        sizeof(int),           &nb,                                VALUE,
        PARSEC_DTD_ARG_END );
#else
    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zunmqr_parsec, options->priority, "unmqr",
        sizeof(cham_side_t),    &side,                              VALUE,
        sizeof(cham_trans_t),    &trans,                             VALUE,
        sizeof(int),           &m,                                 VALUE,
        sizeof(int),           &n,                                 VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,   NULL,                          SCRATCH,
        sizeof(int),           &nb,                                VALUE,
        PARSEC_DTD_ARG_END );
#endif
}
