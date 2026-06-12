/**
 *
 * @file quark/codelet_zperm_allreduce.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon quark codelets to do the reduction
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

#if defined(CHAMELEON_USE_MPI)

void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved  )
{
    (void)options;
    (void)A;
    (void)Am;
    (void)An;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved  )
{
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_invp_row( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     k,
                                           int                     n )
{
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)k;
    (void)n;
}

void
INSERT_TASK_zperm_allreduce_send_invp_col( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     m,
                                           int                     k )
{
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)m;
    (void)k;
}

void
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
                             cham_dir_t              dir,
                             const CHAM_desc_t      *A,
                             CHAM_desc_t            *U,
                             int                     Um,
                             int                     Un,
                             CHAM_ipiv_t            *ipiv,
                             int                     ipivk,
                             int                     k,
                             int                     n,
                             void                   *ws )
{
    (void)options;
    (void)dir;
    (void)A;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
    (void)ipivk;
    (void)k;
    (void)n;
    (void)ws;
}

#endif /* if defined(CHAMELEON_USE_MPI) */
