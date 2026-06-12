/**
 *
 * @file parsec/codelet_zpivot_allreduce.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon parsec codelets to do the reduction
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"

#if defined(CHAMELEON_USE_MPI)

void INSERT_TASK_zpivot_allreduce( const RUNTIME_option_t *options,
                                   CHAM_desc_t            *A,
                                   CHAM_desc_pivot_t      *pivot,
                                   int                     k,
                                   int                     h,
                                   int                     n,
                                   void                   *ws )
{
    (void)options;
    (void)A;
    (void)pivot;
    (void)k;
    (void)h;
    (void)n;
    (void)ws;
}

#endif
