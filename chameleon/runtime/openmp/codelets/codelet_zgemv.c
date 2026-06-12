/**
 *
 * @file openmp/codelet_zgemv.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemv OpenMP codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void
INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                   cham_trans_t trans, int m, int n,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *X, int Xm, int Xn, int incX,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileX = X->get_blktile( X, Xm, Xn );
    CHAM_tile_t *tileY = Y->get_blktile( Y, Ym, Yn );

    if ( alpha == 0. ) {
#pragma omp task firstprivate( m, n, alpha, tileY ) depend( inout:tileY[0] )
        TCORE_zlascal( ChamUpperLower, m, n, alpha, tileA );
    }
    else if ( beta == 0. ) {
#pragma omp task firstprivate( trans, m, n, alpha, tileA, tileX, incX, beta, tileY, incY ) depend( in:tileA[0], tileX[0] ) depend( out:tileY[0] )
        TCORE_zgemv( trans,  m, n,
                     alpha, tileA, tileX, incX,
                     beta,  tileY, incY );
    }
    else {
#pragma omp task firstprivate( trans, m, n, alpha, tileA, tileX, incX, beta, tileY, incY ) depend( in:tileA[0], tileX[0] ) depend( inout:tileY[0] )
        TCORE_zgemv( trans,  m, n,
                     alpha, tileA, tileX, incX,
                     beta,  tileY, incY );
    }
    (void)options;
}
