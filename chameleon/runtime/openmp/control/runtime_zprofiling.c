/**
 *
 * @file openmp/runtime_zprofiling.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP CHAMELEON_Complex64_t kernel progiling
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#include "chameleon_openmp.h"

void RUNTIME_zdisplay_allprofile()
{
    chameleon_warning("RUNTIME_zdisplay_allprofile(OpenMP)", "Profiling is not available with OpenMP");
}

void RUNTIME_zdisplay_oneprofile( cham_tasktype_t kernel )
{
    (void)kernel;
    chameleon_warning("RUNTIME_zdisplay_oneprofile(OpenMP)", "Profiling is not available with OpenMP\n");
}

