/**
 *
 * @file cuda_zgessm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgessm GPU kernel
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
int CUDA_zgessm(
        char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib,
        magma_int_t *ipiv,
        cuDoubleComplex *dL1, magma_int_t lddl1,
        cuDoubleComplex *dL,  magma_int_t lddl,
        cuDoubleComplex *dA,  magma_int_t ldda,
        magma_int_t *info)
{
    int ret;
    /* The kernel is just using the inverted part or nothing */
    if ( lddl1 >= 2*ib ) {
      dL1 += ib;
      ret = magma_zgessm_gpu( storev, m, n, k, ib,
                  ipiv, dL1, lddl1, dL, lddl, dA, ldda, info );
    }
    else {
      ret = magma_zgessm_gpu( storev, m, n, k, ib,
              ipiv, NULL, 1, dL, lddl, dA, ldda, info );
    }

    if (ret != MAGMA_SUCCESS) {
        fprintf(stderr, "Error in MAGMA: %d\n", ret);
        exit(-1);
    }

    return CHAMELEON_SUCCESS;
}
#endif
#endif
