/**
 *
 * @file openmp/codelet_zplssq.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq OpenMP codelet
 *
 * @version 1.4.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *IN,  int INm,  int INn,
                         const CHAM_desc_t *OUT, int OUTm, int OUTn )
{
    CHAM_tile_t *tileIN  = IN->get_blktile( IN, INm, INn );
    CHAM_tile_t *tileOUT = OUT->get_blktile( OUT, OUTm, OUTn );

    assert( tileIN->format  & CHAMELEON_TILE_FULLRANK );
    assert( tileOUT->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( storev, M, N ) depend( in: tileIN[0] ) depend( inout: tileOUT[0] )
    CORE_zplssq( storev, M, N, tileIN->mat, tileOUT->mat );

    (void)options;
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    CHAM_tile_t *tileRESULT = RESULT->get_blktile( RESULT, RESULTm, RESULTn );

    assert( tileRESULT->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( N ) depend( inout: tileRESULT[0] )
    CORE_zplssq2( N, tileRESULT->mat );

    (void)options;
}
