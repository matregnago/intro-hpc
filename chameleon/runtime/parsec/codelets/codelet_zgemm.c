/**
 *
 * @file parsec/codelet_zgemm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm PaRSEC codelet
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Florent Pruvost
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
CORE_zgemm_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int ldc;

    parsec_dtd_unpack_args(
        this_task, &transA, &transB, &m, &n, &k, &alpha, &A, &lda, &B, &ldb, &beta, &C, &ldc );

    CORE_zgemm( transA, transB, m, n, k,
                alpha, A, lda,
                       B, ldb,
                beta,  C, ldc );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
/* One cuBLAS handle per worker thread; persists across tasks (handle creation
 * is non-trivial), but the stream is rebound per task to the GPU exec stream
 * supplied by PaRSEC. Handles leak on thread exit — acceptable for now. */
static __thread cublasHandle_t cham_parsec_cublas_handle = NULL;

static int
CORE_zgemm_parsec_cuda( parsec_execution_stream_t *es,
                        parsec_task_t             *this_task,
                        void                      *cuda_stream )
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m, n, k, lda, ldb, ldc;
    CHAMELEON_Complex64_t alpha, beta;
    CHAMELEON_Complex64_t *A, *B, *C;

    parsec_dtd_unpack_args(
        this_task, &transA, &transB, &m, &n, &k, &alpha, &A, &lda, &B, &ldb, &beta, &C, &ldc );

    /* parsec_dtd_unpack_args returns this_task->data[i].data_in (the host
     * copy). For GPU dispatch, the PaRSEC GPU scheduler has staged each tile
     * to the device and stored the device-side parsec_data_copy_t in
     * this_task->data[i].data_out. Overwrite A/B/C with the device pointers
     * so cuBLAS receives valid GPU addresses. */
    A = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[0].data_out);
    B = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[1].data_out);
    C = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[2].data_out);

    if( NULL == cham_parsec_cublas_handle ) {
        cublasStatus_t cs = cublasCreate( &cham_parsec_cublas_handle );
        if( cs != CUBLAS_STATUS_SUCCESS ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
    }
    cublasSetStream( cham_parsec_cublas_handle, (cudaStream_t)cuda_stream );

    CUDA_zgemm( transA, transB, m, n, k,
                (const cuDoubleComplex *)&alpha,
                (const cuDoubleComplex *)A, lda,
                (const cuDoubleComplex *)B, ldb,
                (const cuDoubleComplex *)&beta,
                (cuDoubleComplex *)C, ldc,
                cham_parsec_cublas_handle );

    (void)es;
    return PARSEC_HOOK_RETURN_DONE;
}
#endif /* CHAMELEON_USE_CUDA && !CHAMELEON_SIMULATION */

void
INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                   cham_trans_t transA, cham_trans_t transB,
                   int m, int n, int k, int nb,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *B, int Bm, int Bn,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    parsec_dtd_taskpool_insert_task_cuda(
        PARSEC_dtd_taskpool,
        CORE_zgemm_parsec, CORE_zgemm_parsec_cuda,
        options->priority, "gemm",
        sizeof(cham_trans_t),    &transA,                           VALUE,
        sizeof(cham_trans_t),    &transB,                           VALUE,
        sizeof(int),           &m,                                VALUE,
        sizeof(int),           &n,                                VALUE,
        sizeof(int),           &k,                                VALUE,
        sizeof(CHAMELEON_Complex64_t),           &alpha,              VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INPUT,
        sizeof(int), &(tileB->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t),           &beta,               VALUE,
        PASSED_BY_REF,     RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        PARSEC_DTD_ARG_END );
#else
    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgemm_parsec, options->priority, "gemm",
        sizeof(cham_trans_t),    &transA,                           VALUE,
        sizeof(cham_trans_t),    &transB,                           VALUE,
        sizeof(int),           &m,                                VALUE,
        sizeof(int),           &n,                                VALUE,
        sizeof(int),           &k,                                VALUE,
        sizeof(CHAMELEON_Complex64_t),           &alpha,              VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INPUT,
        sizeof(int), &(tileB->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t),           &beta,               VALUE,
        PASSED_BY_REF,     RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | INOUT | AFFINITY,
        sizeof(int), &(tileC->ld), VALUE,
        PARSEC_DTD_ARG_END );
#endif

    (void)nb;
}

void
INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                      const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    INSERT_TASK_zgemm( options, transA, transB, m, n, k, nb,
                       alpha, A, Am, An, B, Bm, Bn,
                       beta, C, Cm, Cn );
}
