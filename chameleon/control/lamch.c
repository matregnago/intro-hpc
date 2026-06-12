/**
 *
 * @file lamch.c
 *
 * @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon lamch wrapper to factorize the test with simulation compilations
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-09-18
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/lapacke.h"
#endif

/**
 *
 * Wrapper to the LAPACKE_slamch_work('e') call that returns a pre-defined value
 * if compile with SimGrid.
 *
 * @return The epsilon value for single precision
 *
 */
float
CHAMELEON_slamch( void )
{
#if !defined(CHAMELEON_SIMULATION)
    return LAPACKE_slamch_work( 'e' );
#else
    return 1.e-7;
#endif
}

/**
 *
 * Wrapper to the LAPACKE_dlamch_work('e') call that returns a pre-defined value
 * if compile with SimGrid.
 *
 * @return The epsilon value for single precision
 *
 */
double
CHAMELEON_dlamch( void )
{
#if !defined(CHAMELEON_SIMULATION)
    return LAPACKE_dlamch_work( 'e' );
#else
    return 1.e-15;
#endif
}
