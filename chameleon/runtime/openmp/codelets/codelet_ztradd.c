/**
 *
 * @file openmp/codelet_ztradd.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    if ( alpha == 0. ) {
#pragma omp task firstprivate( uplo, m, n, beta, tileB ) depend( inout:tileB[0] )
        TCORE_zlascal( uplo, m, n, beta, tileB );
    }
    else if ( beta == 0. ) {
#pragma omp task firstprivate( uplo, trans, m, n, alpha, tileA, tileB ) depend( in:tileA[0] ) depend( out:tileB[0] )
        TCORE_ztradd( uplo, trans, m, n, alpha, tileA, beta, tileB );
    }
    else {
#pragma omp task firstprivate( uplo, trans, m, n, alpha, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
        TCORE_ztradd( uplo, trans, m, n, alpha, tileA, beta, tileB );
    }

    (void)options;
    (void)nb;
}
