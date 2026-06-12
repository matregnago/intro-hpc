/**
 *
 * @file parsec/codelet_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpmqrt PaRSEC codelet
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

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "gpucublas.h"
#endif

static inline int
CORE_ztpmqrt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    const CHAMELEON_Complex64_t *V;
    int ldv;
    const CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t *WORK;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &M, &N, &K, &L, &ib, &V, &ldv, &T, &ldt, &A, &lda, &B, &ldb, &WORK );

    CORE_ztpmqrt( side, trans, M, N, K, L, ib,
                  V, ldv, T, ldt, A, lda, B, ldb, WORK );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
/* Per-thread cuBLAS handle + GPU scratch. The scratch is 3*ib*nb elements
 * (same sizing StarPU's cl_ztpmqrt_cuda_func documents). Tile size is fixed
 * across a run, so the realloc fires once per worker. Same lifecycle
 * compromise as codelet_zpotrf.c: handle and scratch leak at thread exit. */
static __thread cublasHandle_t   cham_parsec_tpmqrt_cublas_handle = NULL;
static __thread cuDoubleComplex *cham_parsec_tpmqrt_scratch = NULL;
static __thread int              cham_parsec_tpmqrt_scratch_elts = 0;

static int
CORE_ztpmqrt_parsec_cuda( parsec_execution_stream_t *es,
                          parsec_task_t             *this_task,
                          void                      *cuda_stream )
{
    cham_side_t side;
    cham_trans_t trans;
    int M, N, K, L, ib;
    const CHAMELEON_Complex64_t *V, *T;
    int ldv, ldt, lda, ldb;
    CHAMELEON_Complex64_t *A, *B, *WORK_unused;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &M, &N, &K, &L, &ib,
        &V, &ldv, &T, &ldt, &A, &lda, &B, &ldb, &WORK_unused );

    V = (const CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[0].data_out);
    T = (const CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[1].data_out);
    A = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[2].data_out);
    B = (CHAMELEON_Complex64_t *)parsec_data_copy_get_ptr(this_task->data[3].data_out);

    if( NULL == cham_parsec_tpmqrt_cublas_handle ) {
        cublasStatus_t cs = cublasCreate( &cham_parsec_tpmqrt_cublas_handle );
        if( cs != CUBLAS_STATUS_SUCCESS ) {
            return PARSEC_HOOK_RETURN_ERROR;
        }
    }
    cublasSetStream( cham_parsec_tpmqrt_cublas_handle, (cudaStream_t)cuda_stream );

    /* CUDA_ztpmqrt needs 3*ib*nb workspace; nb is implicit in the tile
     * dimension (N when side==Left, M when side==Right). Use max(M,N) as an
     * upper bound that covers both cases without parsing side. */
    int nb_upper = (M > N) ? M : N;
    int need = 3 * ib * nb_upper;
    if( need > cham_parsec_tpmqrt_scratch_elts ) {
        if( cham_parsec_tpmqrt_scratch != NULL ) {
            cudaFree( cham_parsec_tpmqrt_scratch );
        }
        if( cudaMalloc( (void **)&cham_parsec_tpmqrt_scratch,
                        (size_t)need * sizeof(cuDoubleComplex) ) != cudaSuccess ) {
            cham_parsec_tpmqrt_scratch = NULL;
            cham_parsec_tpmqrt_scratch_elts = 0;
            return PARSEC_HOOK_RETURN_ERROR;
        }
        cham_parsec_tpmqrt_scratch_elts = need;
    }

    CUDA_ztpmqrt( side, trans, M, N, K, L, ib,
                  (const cuDoubleComplex *)V, ldv,
                  (const cuDoubleComplex *)T, ldt,
                  (cuDoubleComplex *)A, lda,
                  (cuDoubleComplex *)B, ldb,
                  cham_parsec_tpmqrt_scratch, need,
                  cham_parsec_tpmqrt_cublas_handle );

    (void)es;
    (void)WORK_unused;
    return PARSEC_HOOK_RETURN_DONE;
}
#endif /* CHAMELEON_USE_CUDA && !CHAMELEON_SIMULATION */

void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const CHAM_desc_t *V, int Vm, int Vn,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    parsec_dtd_taskpool_insert_task_cuda(
        PARSEC_dtd_taskpool,
        CORE_ztpmqrt_parsec, CORE_ztpmqrt_parsec_cuda,
        options->priority, "tpmqrt",
        sizeof(cham_side_t), &side,  VALUE,
        sizeof(cham_trans_t), &trans, VALUE,
        sizeof(int),        &M,     VALUE,
        sizeof(int),        &N,     VALUE,
        sizeof(int),        &K,     VALUE,
        sizeof(int),        &L,     VALUE,
        sizeof(int),        &ib,    VALUE,
        PASSED_BY_REF,       RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), chameleon_parsec_get_arena_index( V ) | INPUT,
        sizeof(int), &(tileV->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL, SCRATCH,
        PARSEC_DTD_ARG_END );
#else
    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztpmqrt_parsec, options->priority, "tpmqrt",
        sizeof(cham_side_t), &side,  VALUE,
        sizeof(cham_trans_t), &trans, VALUE,
        sizeof(int),        &M,     VALUE,
        sizeof(int),        &N,     VALUE,
        sizeof(int),        &K,     VALUE,
        sizeof(int),        &L,     VALUE,
        sizeof(int),        &ib,    VALUE,
        PASSED_BY_REF,       RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), chameleon_parsec_get_arena_index( V ) | INPUT,
        sizeof(int), &(tileV->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | INPUT,
        sizeof(int), &(tileT->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INOUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL, SCRATCH,
        PARSEC_DTD_ARG_END );
#endif
}
