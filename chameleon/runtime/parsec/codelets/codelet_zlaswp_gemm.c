/**
 *
 * @file parsec/codelet_zlaswp_gemm.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon parsec codelets to apply zlaswp_set and gemm on a panel
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswp_gemm( const RUNTIME_option_t *options,
                              cham_side_t             side,
                              cham_dir_t              dir,
                              int                     m0,
                              int                     m,
                              int                     n,
                              int                     k,
                              void                   *ws,
                              const CHAM_ipiv_t      *ipiv, int ipivk,
                              const CHAM_desc_t      *WA,   int WAm, int WAn,
                              const CHAM_desc_t      *A,    int Am,  int An,
                              const CHAM_desc_t      *B,    int Bm,  int Bn,
                              const CHAM_desc_t      *C,    int Cm,  int Cn,
                              void                  **clargs_ptr )
{
    (void)options;
    (void)side;
    (void)dir;
    (void)m0;
    (void)m;
    (void)n;
    (void)k;
    (void)ws;
    (void)ipiv;
    (void)ipivk;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)A;
    (void)Am;
    (void)An;
    (void)B;
    (void)Bm;
    (void)Bn;
    (void)C;
    (void)Cm;
    (void)Cn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_gemm_flush( const RUNTIME_option_t *options,
                                    cham_dir_t              dir,
                                    const CHAM_ipiv_t      *ipiv, int ipivk,
                                    const CHAM_desc_t      *WA,   int WAm, int WAn,
                                    const CHAM_desc_t      *B,    int Bm,  int Bn,
                                    void                  **clargs_ptr )
{
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)B;
    (void)Bm;
    (void)Bn;
    (void)clargs_ptr;
}
