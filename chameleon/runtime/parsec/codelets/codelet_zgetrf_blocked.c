/**
 *
 * @file parsec/codelet_zgetrf_blocked.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_blocked Parsec codelets
 *
 * @version 1.4.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib, int readUp,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t       *ipiv,
                                      CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)ib;
    (void)readUp;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
    (void)pivot;
}

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib, int readUp,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)ib;
    (void)readUp;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)pivot;
}
