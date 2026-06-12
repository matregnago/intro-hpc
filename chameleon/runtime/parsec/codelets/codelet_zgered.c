/**
 *
 * @file parsec/codelet_zgered.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgered PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Ana Hourcau
 * @date 2024-07-17
 * @precisions normal z -> d
 *
 */
#include "chameleon_parsec.h"

void INSERT_TASK_zgered( const RUNTIME_option_t *options,
                         double threshold, int m, int n,
                         const CHAM_desc_t *A,     int Am,  int An,
                         const CHAM_desc_t *Wnorm, int Wnm, int Wnn )
{
    fprintf( stderr, "WARNING: gered kernel is not available with PaRSEC\n" );

    (void)options;
    (void)threshold;
    (void)m;
    (void)n;
    (void)A;
    (void)Am;
    (void)An;
    (void)Wnorm;
    (void)Wnm;
    (void)Wnn;
}
