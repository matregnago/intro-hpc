/**
 *
 * @file quark/codelet_zperm_reduce.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon quark codelets to do the reduction
 *
 * @version 1.4.0
 * @author Matteo Marcox
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

void
INSERT_TASK_zperm_reduce( const RUNTIME_option_t *options,
                          cham_dir_t              dir,
                          const CHAM_desc_t      *A,
                          int                     m,
                          int                     n,
                          CHAM_ipiv_t            *ipiv,
                          int                     ipivk,
                          const CHAM_desc_t      *Wu,
                          int                     Wum,
                          int                     Wun,
                          void                   *ws,
                          int                     Wm,
                          int                     Wn )
{
    (void)options;
    (void)dir;
    (void)A;
    (void)m;
    (void)n;
    (void)ipiv;
    (void)ipivk;
    (void)Wu;
    (void)Wum;
    (void)Wun;
    (void)ws;
    (void)Wm;
    (void)Wn;
}
