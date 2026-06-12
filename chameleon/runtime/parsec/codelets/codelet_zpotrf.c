/**
 *
 * @file parsec/codelet_zpotrf.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf PaRSEC codelet
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
#include <cusolverDn.h>
#include "gpucublas.h"
#endif

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
static inline int
CORE_zpotrf_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    int tempkm, ldak, iinfo, info;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    CHAMELEON_Complex64_t *A;

    parsec_dtd_unpack_args(
        this_task, &uplo, &tempkm, &A, &ldak, &iinfo, &sequence, &request );

    CORE_zpotrf( uplo, tempkm, A, ldak, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    (void)info;
    (void)iinfo;
    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
/* Per-thread cuSOLVER state. Handle, workspace, and devInfo are created
 * lazily on first dispatch and reused for the lifetime of the worker thread.
 * Same lifecycle compromise as the cuBLAS handle in codelet_zgemm.c: handles
 * and device allocations leak at thread exit. */
static __thread cusolverDnHandle_t cham_parsec_potrf_cusolver_handle = NULL;
static __thread cuDoubleComplex   *cham_parsec_potrf_workspace = NULL;
static __thread int                cham_parsec_potrf_workspace_elts = 0;
static __thread int               *cham_parsec_potrf_d_info = NULL;

static int
CORE_zpotrf_parsec_cuda( parsec_execution_stream_t *es,
                         parsec_task_t             *this_task,
                         void                      *cuda_stream )
{
    cham_uplo_t uplo;
    int tempkm, ldak, iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    CHAMELEON_Complex64_t *A;

    parsec_dtd_unpack_args(
        this_task, &uplo, &tempkm, &A, &ldak, &iinfo, &sequence, &request );

    A = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[0].data_out);

    cudaStream_t stream = (cudaStream_t)cuda_stream;

    if( NULL == cham_parsec_potrf_cusolver_handle ) {
        cusolverStatus_t cs = cusolverDnCreate( &cham_parsec_potrf_cusolver_handle );
        if( cs != CUSOLVER_STATUS_SUCCESS ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
        if( cudaMalloc( (void **)&cham_parsec_potrf_d_info, sizeof(int) ) != cudaSuccess ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
    }
    cusolverDnSetStream( cham_parsec_potrf_cusolver_handle, stream );

    /* Size workspace from cuSOLVER; the tile size is fixed across the run so
     * this reallocates exactly once. */
    int lwork = 0;
    if( cusolverDnZpotrf_bufferSize( cham_parsec_potrf_cusolver_handle,
                                     chameleon_cublas_const(uplo),
                                     tempkm,
                                     (cuDoubleComplex *)A, ldak,
                                     &lwork ) != CUSOLVER_STATUS_SUCCESS ) {
        return PARSEC_HOOK_RETURN_ERROR;
    }
    if( lwork > cham_parsec_potrf_workspace_elts ) {
        if( cham_parsec_potrf_workspace != NULL ) {
            cudaFree( cham_parsec_potrf_workspace );
        }
        if( cudaMalloc( (void **)&cham_parsec_potrf_workspace,
                        (size_t)lwork * sizeof(cuDoubleComplex) ) != cudaSuccess ) {
            cham_parsec_potrf_workspace = NULL;
            cham_parsec_potrf_workspace_elts = 0;
            return PARSEC_HOOK_RETURN_ERROR;
        }
        cham_parsec_potrf_workspace_elts = lwork;
    }

    CUDA_zpotrf( uplo, tempkm,
                 (cuDoubleComplex *)A, ldak,
                 cham_parsec_potrf_workspace, lwork,
                 cham_parsec_potrf_d_info,
                 cham_parsec_potrf_cusolver_handle );

    /* Mirror StarPU's codelet_zpotrf cuda body: pull devInfo back so that
     * non-positive-definite tiles surface through RUNTIME_sequence_flush.
     * Forces a stream sync — unavoidable for this kernel. */
    int info = 0;
    cudaMemcpyAsync( &info, cham_parsec_potrf_d_info, sizeof(int),
                     cudaMemcpyDeviceToHost, stream );
    cudaStreamSynchronize( stream );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)es;
    (void)iinfo;
    return PARSEC_HOOK_RETURN_DONE;
}
#endif /* CHAMELEON_USE_CUDA && !CHAMELEON_SIMULATION */

void INSERT_TASK_zpotrf(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    parsec_dtd_taskpool_insert_task_cuda(
        PARSEC_dtd_taskpool,
        CORE_zpotrf_parsec, CORE_zpotrf_parsec_cuda,
        options->priority, "potrf",
        sizeof(cham_uplo_t),                 &uplo,                             VALUE,
        sizeof(int),                 &n,                                VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),                 &iinfo,                            VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),              VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),               VALUE,
        PARSEC_DTD_ARG_END );
#else
    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zpotrf_parsec, options->priority, "potrf",
        sizeof(cham_uplo_t),                 &uplo,                             VALUE,
        sizeof(int),                 &n,                                VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),                 &iinfo,                            VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),              VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),               VALUE,
        PARSEC_DTD_ARG_END );
#endif

    (void)nb;
}
