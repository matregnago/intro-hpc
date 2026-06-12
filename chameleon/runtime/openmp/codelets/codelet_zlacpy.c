/**
 *
 * @file openmp/codelet_zlacpy.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( uplo, m, n, tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        TCORE_zlacpy( uplo, m, n, tileA, tileB );
    }

    (void)options;
}

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( uplo, m, n, displA, tileA, lda, displB, tileB, ldb ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        TCORE_zlacpyx( uplo, m, n, displA, tileA, lda, displB, tileB, ldb );
    }

    (void)options;
}
