/**
 *
 * @file core_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlauum CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Computes the product U * U^H or L^H * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If uplo = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If uplo = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.

 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 *
 * @param[in] N
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L.
 *          On exit, if uplo = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U^H;
 *          if uplo = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L^H * L.

 *
 * @param[in] LDA
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 * @param[out] info
 *          - 0 on successful exit
 *          - < 0 if -i, the i-th argument had an illegal value
 *
 */
void CORE_zlauum(cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, int LDA)
{
    LAPACKE_zlauum_work(LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
}


