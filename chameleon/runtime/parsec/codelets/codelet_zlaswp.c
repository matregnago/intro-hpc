/**
 *
 * @file parsec/codelet_zlaswp.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC codelets to apply zlaswp on a panel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlaswp_get_parsec( parsec_execution_stream_t *context,
                        parsec_task_t             *this_task )
{
    cham_side_t            side;
    int                    m0, m, n, k, lda, ldb, *perm;
    CHAMELEON_Complex64_t *A, *B;

    parsec_dtd_unpack_args( this_task, &side, &m0, &m, &n, &k, &A, &lda, &B, &ldb, &perm );

    CORE_zlaswp_get( side, m0, m, n, k, A, lda, B, ldb, perm );
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn )
{
    assert( 0 );
    (void)options;
    (void)side;
    (void)dir;
    (void)m0;
    (void)m;
    (void)n;
    (void)k;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Am;
    (void)An;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
}

static inline int
CORE_zlaswp_set_parsec( parsec_execution_stream_t *context,
                        parsec_task_t             *this_task )
{
    cham_side_t            side;
    int                    m0, m, n, k, lda, ldb, *invp;
    CHAMELEON_Complex64_t *A, *B;

    parsec_dtd_unpack_args( this_task, &side, &m0, &m, &n, &k, &A, &lda, &B, &ldb, &invp );

    CORE_zlaswp_set( side, m0, m, n, k, A, lda, B, ldb, invp );
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t              dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *WA, int WAm, int WAn,
                             const CHAM_desc_t *A,  int Am,  int An )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileWA = WA->get_blktile( WA, WAm, WAn );
    CHAM_tile_t *tileA  = A->get_blktile(  A,  Am,  An  );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlaswp_set_parsec, options->priority, "laswp_set",
        sizeof(cham_side_t),  &side,   VALUE,
        sizeof(int),          &m0,     VALUE,
        sizeof(int),          &m,      VALUE,
        sizeof(int),          &n,      VALUE,
        sizeof(int),          &k,      VALUE,
        PASSED_BY_REF, RTBLKADDR(WA, ChamComplexDouble, WAm, WAn), chameleon_parsec_get_arena_index( WA )        | INPUT,
        sizeof(int),         &(tileWA->ld), VALUE,
        PASSED_BY_REF, RTBLKADDR(A,  ChamComplexDouble, Am,  An ), chameleon_parsec_get_arena_index( A )         | INOUT,
        sizeof(int),         &(tileA->ld), VALUE,
        PASSED_BY_REF, RUNTIME_ipiv_getinvp( ipiv, ipivk ),        chameleon_parsec_get_arena_index_invp( ipiv ) | INPUT,
        PARSEC_DTD_ARG_END );

    (void)dir;
}

#if defined(CHAMELEON_USE_MPI)

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An )
{
    assert( 0 );
    (void)options;
    (void)ws;
    (void)Wm;
    (void)Wn;
    (void)A;
    (void)Am;
    (void)An;
}

#endif /* if defined(CHAMELEON_USE_MPI) */
