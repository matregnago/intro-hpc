/**
 *
 * @file core_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgetrf CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

int CORE_zgetrf( int m, int n,
                 CHAMELEON_Complex64_t *A, int lda,
                 int *IPIV, int *info )
{
    *info = LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    return CHAMELEON_SUCCESS;
}
