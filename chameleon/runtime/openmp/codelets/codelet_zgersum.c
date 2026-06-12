/**
 *
 * @file openmp/codelet_zgersum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgersum OpenMP codelet
 *
 * @version 1.4.0
 * @author Romain Peressoni
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"

void
RUNTIME_zgersum_set_methods( const CHAM_desc_t *A, int Am, int An )
{
    fprintf( stderr, "WARNING: Reductions are not available with OpenMP\n" );

    (void)A;
    (void)Am;
    (void)An;
}

void
RUNTIME_zgersum_submit_tree( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An )
{
    (void)options;
    (void)A;
    (void)Am;
    (void)An;
}
