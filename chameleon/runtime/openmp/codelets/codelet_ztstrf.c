/**
 *
 * @file openmp/codelet_ztstrf.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *U, int Um, int Un,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *L, int Lm, int Ln,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileU = U->get_blktile( U, Um, Un );
    CHAM_tile_t *tileL = L->get_blktile( L, Lm, Ln );
    int ws_size = options->ws_wsize;
#pragma omp task firstprivate( ws_size, m, n, ib, nb, tileU, tileA, tileL, IPIV, iinfo ) depend( inout:tileA[0], tileU[0], tileL[0] )
    {
      CHAMELEON_Complex64_t work[ws_size];
      TCORE_ztstrf( m, n, ib, nb, tileU, tileA, tileL, IPIV, work, nb, &iinfo );
    }

    (void)options;
    (void)nb;
    (void)check_info;
    (void)iinfo;
}
