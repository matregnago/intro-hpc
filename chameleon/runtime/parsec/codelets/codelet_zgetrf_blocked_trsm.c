/**
 *
 * @file starpu/codelet_zgetrf_blocked_trsm.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel PaRSEC codelets
 *
 * @version 1.4.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zgetrf_cpy_pivrow_in_Up( const RUNTIME_option_t *options,
                                          CHAM_desc_t            *Up,
                                          int                     Upm,
                                          int                     k,
                                          int                     h,
                                          int                     ib,
                                          int                     n,
                                          CHAM_desc_pivot_t      *pivot )
{
    (void)options;
    (void)Up;
    (void)Upm;
    (void)k;
    (void)h;
    (void)ib;
    (void)n;
    (void)pivot;
}

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int                     m,
                                      int                     n,
                                      int                     h,
                                      int                     ib,
                                      CHAM_desc_t            *U,
                                      int                     Um,
                                      int                     Un,
                                      CHAM_desc_pivot_t      *pivot )
{
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)ib;
    (void)U;
    (void)Um;
    (void)Un;
    (void)pivot;
}
