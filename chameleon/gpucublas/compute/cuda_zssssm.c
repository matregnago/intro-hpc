/**
 *
 * @file cuda_zssssm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zssssm GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#if defined(CHAMELEON_USE_MAGMA)
#if defined(HAVE_MAGMA_GETRF_INCPIV_GPU)
int CUDA_zssssm(
        magma_storev_t storev, magma_int_t m1, magma_int_t n1,
        magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib,
        magmaDoubleComplex *dA1, magma_int_t ldda1,
        magmaDoubleComplex *dA2, magma_int_t ldda2,
        magmaDoubleComplex *dL1, magma_int_t lddl1,
        magmaDoubleComplex *dL2, magma_int_t lddl2,
        magma_int_t *IPIV, magma_int_t *info)
{
    magma_zssssm_gpu(
        storev, m1, n1, m2, n2, k, ib,
        dA1, ldda1, dA2, ldda2,
        dL1, lddl1, dL2, lddl2,
        IPIV, info);
    return CHAMELEON_SUCCESS;
}
#endif
#endif
