/**
 *
 * @file openmp/codelet_zgemm.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm OpenMP codelet
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

void
INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                   cham_trans_t transA, cham_trans_t transB,
                   int m, int n, int k, int nb,
                   CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                const CHAM_desc_t *B, int Bm, int Bn,
                   CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

    if ( alpha == 0. ) {
#pragma omp task firstprivate( m, n, beta, tileC ) depend( inout:tileC[0] )
        TCORE_zlascal( ChamUpperLower, m, n, beta, tileC );
    }
    else if ( beta == 0. ) {
#pragma omp task firstprivate( transA, transB, m, n, k, alpha, tileA, tileB, beta, tileC ) depend( in:tileA[0], tileB[0] ) depend( out:tileC[0] )
        TCORE_zgemm( transA, transB,
                     m, n, k,
                     alpha, tileA,
                     tileB,
                     beta, tileC );
    }
    else {
#pragma omp task firstprivate( transA, transB, m, n, k, alpha, tileA, tileB, beta, tileC ) depend( in:tileA[0], tileB[0] ) depend( inout:tileC[0] )
        TCORE_zgemm( transA, transB,
                     m, n, k,
                     alpha, tileA,
                     tileB,
                     beta, tileC );
    }

    (void)options;
    (void)nb;
}

void
INSERT_TASK_zgemm_Astat( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    INSERT_TASK_zgemm( options, transA, transB, m, n, k, nb,
                       alpha, A, Am, An, B, Bm, Bn,
                       beta, C, Cm, Cn );
}
