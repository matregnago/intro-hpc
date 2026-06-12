/**
 *
 * @file cuda_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztstrf GPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

#if defined(CHAMELEON_USE_MAGMA) && 0
int CUDA_ztstrf(
        char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
        cuDoubleComplex *hU, magma_int_t ldhu, cuDoubleComplex *dU, magma_int_t lddu,
        cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
        cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
        magma_int_t *ipiv,
        cuDoubleComplex *hwork, magma_int_t ldhwork,
        cuDoubleComplex *dwork, magma_int_t lddwork,
        magma_int_t *info)
{

    magma_ztstrf_gpu( storev, m, n, ib, nb,
                      hU, ldhu, dU, lddu,
                      hA, ldha, dA, ldda,
                      hL, ldhl, dL, lddl,
                      ipiv,
                      hwork, ldhwork, dwork, lddwork,
                      info );
    return CHAMELEON_SUCCESS;
}
#endif
