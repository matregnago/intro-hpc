/**
 *
 * @file starpu/codelet_zgetrf_batched.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel batched OpenMP codelets
 *
 * @version 1.4.0
 * @comment Codelets to perform batched panel factorization with partial pivoting
 *
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void
INSERT_TASK_zgetrf_panel_offdiag_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          void **clargs_ptr,
                                          CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)A;
    (void)Am;
    (void)An;
    (void)clargs_ptr;
    (void)pivot;
}

void
INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                void **clargs_ptr,
                                                CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)A;
    (void)An;
    (void)clargs_ptr;
    (void)pivot;
}

void
INSERT_TASK_zgetrf_panel_blocked_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0, int readUp,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          CHAM_desc_t *U, int Um, int Un,
                                          void **clargs_ptr,
                                          CHAM_ipiv_t *ipiv,
                                          CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)readUp;
    (void)ws;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
    (void)ipiv;
    (void)pivot;
}

void
INSERT_TASK_zgetrf_panel_blocked_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An, int readUp,
                                                CHAM_desc_t *U, int Um, int Un,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv,
                                                CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)A;
    (void)An;
    (void)readUp;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
    (void)ipiv;
    (void)pivot;
}
