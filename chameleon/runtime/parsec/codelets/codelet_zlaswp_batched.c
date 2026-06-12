/**
 *
 * @file parsec/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Parsec codelets to apply zlaswp on a panel
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 cham_side_t             side,
                                 cham_dir_t              dir,
                                 int                     m0,
                                 int                     m,
                                 int                     n,
                                 int                     k,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv, int ipivk,
                                 const CHAM_desc_t      *A,   int Am,   int An,
                                 const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                 const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                 void                  **clargs_ptr )
{
    assert( 0 );
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
    (void)A;
    (void)Am;
    (void)An;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                       const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                       void                  **clargs_ptr )
{
    assert( 0 );
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_get_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *A,   int Am,   int An,
                                     const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                     void                  **clargs_ptr )
{
    assert( 0 );
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
    (void)A;
    (void)Am;
    (void)An;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_get_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WAP, int WAPm, int WAPn,
                                           void                  **clargs_ptr )
{
    assert( 0 );
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)WAP;
    (void)WAPm;
    (void)WAPn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_set_batched( const RUNTIME_option_t *options,
                                     cham_side_t             side,
                                     cham_dir_t              dir,
                                     int                     m0,
                                     int                     m,
                                     int                     n,
                                     int                     k,
                                     void                   *ws,
                                     const CHAM_ipiv_t      *ipiv, int ipivk,
                                     const CHAM_desc_t      *A,   int Am,   int An,
                                     const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                     void                  **clargs_ptr )
{
    assert( 0 );
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
    (void)A;
    (void)Am;
    (void)An;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_set_batched_flush( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           const CHAM_ipiv_t      *ipiv, int ipivk,
                                           const CHAM_desc_t      *WA,  int WAm,  int WAn,
                                           void                  **clargs_ptr )
{
    assert( 0 );
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)ipivk;
    (void)WA;
    (void)WAm;
    (void)WAn;
    (void)clargs_ptr;
}
