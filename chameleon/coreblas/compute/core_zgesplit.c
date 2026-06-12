/**
 *
 * @file core_zgesplit.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgesplit CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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
 ***************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zgesplit splits a matrix A into two parts (Upper/Lower), A keeps its
 *  lower/upper part unchanged and the other part is filled with zeros. Ones
 *  can be optionally set on the diagonal.  The part of A which is erased is
 *  copied in B.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in,out] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zgesplit(cham_side_t side, cham_diag_t diag,
                  int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA,
                  CHAMELEON_Complex64_t *B, int LDB)
{
    cham_uplo_t uplo;

    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(7, "Illegal value of LDB");
        return -7;
    }

    if (side == ChamLeft){
        uplo = ChamUpper;
    } else{
        uplo = ChamLower;
    }

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                        chameleon_lapack_const(ChamUpperLower),
                        M, N, A, LDA, B, LDB);

    LAPACKE_zlaset_work(LAPACK_COL_MAJOR,
                        chameleon_lapack_const(uplo),
                        M, N, 0., 1., A, LDA);

    (void)diag;
    return CHAMELEON_SUCCESS;
}
