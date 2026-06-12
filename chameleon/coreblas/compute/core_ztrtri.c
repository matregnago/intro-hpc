/**
 *
 * @file core_ztrtri.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ztrtri CPU kernel
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
 *  Computes the inverse of an upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = ChamNonUnit: A is non-unit triangular;
 *          = ChamUnit:    A is unit triangular.
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If uplo = 'U', the
 *          leading n-by-n upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If uplo =
 *          'L', the leading n-by-n lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If diag =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 * @param[in] info
 *            < 0 if -i, the i-th argument had an illegal value
 *            > 0 if i, A(i,i) is exactly zero.  The triangular
 *            matrix is singular and its inverse can not be computed.
 *
 */
void CORE_ztrtri(cham_uplo_t uplo, cham_diag_t diag, int N, CHAMELEON_Complex64_t *A, int LDA, int *info)
{
    *info = LAPACKE_ztrtri_work(
        LAPACK_COL_MAJOR,
        chameleon_lapack_const(uplo), chameleon_lapack_const(diag),
        N, A, LDA);
}


