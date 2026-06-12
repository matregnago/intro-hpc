/**
 *
 * @file parsec/codelet_ztrsm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm PaRSEC codelet
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
CORE_ztrsm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t trans;
    cham_diag_t diag;
    int tempmm, nb, ldak, ldam;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *T;
    CHAMELEON_Complex64_t *C;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &trans, &diag, &tempmm, &nb, &alpha, &T, &ldak, &C, &ldam );

    CORE_ztrsm( side, uplo, trans, diag,
                tempmm, nb, alpha,
                T, ldak, C, ldam );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
static __thread cublasHandle_t cham_parsec_trsm_cublas_handle = NULL;

static int
CORE_ztrsm_parsec_cuda( parsec_execution_stream_t *es,
                        parsec_task_t             *this_task,
                        void                      *cuda_stream )
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t trans;
    cham_diag_t diag;
    int tempmm, nb, ldak, ldam;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *T;
    CHAMELEON_Complex64_t *C;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &trans, &diag, &tempmm, &nb, &alpha, &T, &ldak, &C, &ldam );

    /* See codelet_zgemm.c: DTD's unpack returns the host data pointer; the
     * staged device copy lives in data_out. Override T and C with the device
     * addresses before handing them to cuBLAS. */
    T = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[0].data_out);
    C = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[1].data_out);

    if( NULL == cham_parsec_trsm_cublas_handle ) {
        cublasStatus_t cs = cublasCreate( &cham_parsec_trsm_cublas_handle );
        if( cs != CUBLAS_STATUS_SUCCESS ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
    }
    cublasSetStream( cham_parsec_trsm_cublas_handle, (cudaStream_t)cuda_stream );

    CUDA_ztrsm( side, uplo, trans, diag,
                tempmm, nb,
                (const cuDoubleComplex *)&alpha,
                (const cuDoubleComplex *)T, ldak,
                (cuDoubleComplex *)C, ldam,
                cham_parsec_trsm_cublas_handle );

    (void)es;
    return PARSEC_HOOK_RETURN_DONE;
}
#endif /* CHAMELEON_USE_CUDA && !CHAMELEON_SIMULATION */

void INSERT_TASK_ztrsm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    parsec_dtd_taskpool_insert_task_cuda(
        PARSEC_dtd_taskpool,
        CORE_ztrsm_parsec, CORE_ztrsm_parsec_cuda,
        options->priority, "Trsm",
        sizeof(cham_side_t),  &side,                     VALUE,
        sizeof(cham_uplo_t),  &uplo,                     VALUE,
        sizeof(cham_trans_t), &transA,                   VALUE,
        sizeof(cham_diag_t),  &diag,                     VALUE,
        sizeof(int),          &m,                        VALUE,
        sizeof(int),          &n,                        VALUE,
        sizeof(CHAMELEON_Complex64_t),           &alpha,      VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
#else
    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztrsm_parsec, options->priority, "Trsm",
        sizeof(cham_side_t),  &side,                     VALUE,
        sizeof(cham_uplo_t),  &uplo,                     VALUE,
        sizeof(cham_trans_t), &transA,                   VALUE,
        sizeof(cham_diag_t),  &diag,                     VALUE,
        sizeof(int),          &m,                        VALUE,
        sizeof(int),          &n,                        VALUE,
        sizeof(CHAMELEON_Complex64_t),           &alpha,      VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );
#endif

    (void)nb;
}
