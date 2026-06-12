/**
 *
 * @file quark/runtime_zlocality.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark CHAMELEON_Complex64_t kernel locality management
 *
 * @version 1.4.0
 * @author Vijay Joshi
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_quark.h"

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
    chameleon_warning("RUNTIME_zlocality_allrestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_onerestrict( cham_tasktype_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
    chameleon_warning("RUNTIME_zlocality_onerestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_allrestore( )
{
    chameleon_warning("RUNTIME_zlocality_allrestore(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_onerestore( cham_tasktype_t kernel )
{
    (void)kernel;
    chameleon_warning("RUNTIME_zlocality_onerestore(quark)", "Kernel locality cannot be specified with Quark");
}
