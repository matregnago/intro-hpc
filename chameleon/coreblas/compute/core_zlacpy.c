/**
 *
 * @file core_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlacpy CPU kernel
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
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - ChamGeneral: entire A,
 *          - ChamUpper:   upper triangle,
 *          - ChamLower:   lower triangle.
 *
 * @param[in] M
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 */
void CORE_zlacpy(cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB)
{
    LAPACKE_zlacpy_work(
        LAPACK_COL_MAJOR,
        chameleon_lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}



