/**
 *
 * @file core_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlag2c CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Converts m-by-n matrix A from double complex to single complex precision.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *          n >= 0.
 *
 * @param[in] A
 *          The lda-by-n matrix in double complex precision to convert.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          On exit, the converted LDB-by-n matrix in single complex precision.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix As.
 *          ldas >= max(1,m).
 *
 */
void CORE_zlag2c(int m, int n,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex32_t *B, int ldb, int *info)
{
    *info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}
/**/

/**/

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Converts m-by-n matrix A from single complex to double complex precision.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix As.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix As.
 *          n >= 0.
 *
 * @param[in] As
 *          The ldas-by-n matrix in single complex precision to convert.
 *
 * @param[in] ldas
 *          The leading dimension of the matrix As.
 *          ldas >= max(1,m).
 *
 * @param[out] A
 *          On exit, the converted lda-by-n matrix in double complex precision.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *          lda >= max(1,m).
 *
 */
void CORE_clag2z(int m, int n,
                 const CHAMELEON_Complex32_t *A, int lda,
                 CHAMELEON_Complex64_t *B, int ldb)
{
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}



