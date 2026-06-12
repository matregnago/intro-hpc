/**
 *
 * @file quark/codelet_zlaswp.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark codelets to apply zlaswp on a panel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static void CORE_zlaswp_get_quark( Quark *quark )
{
    cham_side_t  side;
    int          m0, k, *perm;
    CHAM_tile_t *A, *B;

    quark_unpack_args_6( quark, side, m0, k, perm, A, B );

    TCORE_zlaswp_get( side, m0, A->m, A->n, k, A, B, perm );
}

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A,   int Am,   int An,
                             const CHAM_desc_t *WAP, int WAPm, int WAPn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASWP;

    QUARK_Insert_Task(
        opt->quark, CORE_zlaswp_get_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_side_t),  &side, VALUE,
        sizeof(int),          &m0,   VALUE,
        sizeof(int),          &k,    VALUE,
        sizeof(int*),         RUNTIME_ipiv_getperm( ipiv, ipivk ),           INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(A,   ChamComplexDouble, Am,   An  ), INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(WAP, ChamComplexDouble, WAPm, WAPn), INOUT,
        0 );

    (void)dir;
    (void)m;
    (void)n;
}

static void CORE_zlaswp_set_quark( Quark *quark )
{
    cham_side_t  side;
    int          m0, m, n, k, *invp;
    CHAM_tile_t *A, *B;

    quark_unpack_args_8( quark, side, m0, m, n, k, invp, A, B );

    TCORE_zlaswp_set( side, m0, A->m, A->n, k, A, B, invp );
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t              dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *WA, int WAm, int WAn,
                             const CHAM_desc_t *A,  int Am,  int An )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASWP;

    QUARK_Insert_Task(
        opt->quark, CORE_zlaswp_set_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_side_t),  &side, VALUE,
        sizeof(int),          &m0,   VALUE,
        sizeof(int),          &m,    VALUE,
        sizeof(int),          &n,    VALUE,
        sizeof(int),          &k,    VALUE,
        sizeof(int*),         RUNTIME_ipiv_getinvp( ipiv, ipivk ),        INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(WA, ChamComplexDouble, WAm, WAn), INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(A,  ChamComplexDouble, Am, An),   INOUT,
        0 );

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
