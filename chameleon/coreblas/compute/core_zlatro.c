/**
 *
 * @file core_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlatro CPU kernel
 *
 * @version 1.4.0
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/lapacke.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zlatro transposes a m-by-n matrix out of place.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = ChamUpper: the upper triangle of A and the lower triangle of B
 *          are referenced.
 *          = ChamLower: the lower triangle of A and the upper triangle of B
 *          are referenced.
 *          = ChamUpperLower: All A and B are referenced.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          = ChamNoTrans:   B is a copy of A (equivalent to zlacpy);
 *          = ChamTrans:     B is the transpose of A;
 *          = ChamConjTrans: B is the conjugate transpose of A.
 *
 * @param[in] M
 *         Number of rows of the matrix A and number of columns of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of rows of the matrix A and the matrix B, if trans == PasmaNoTrans.
 *
 * @param[in] N
 *         Number of columns of the matrix A and number of rows of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of columns of the matrix A and of the matrix B, if trans == ChamNoTrans.
 *
 * @param[in] A
 *         Matrix of size LDA-by-N, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDA-by-M, if trans == PasmaNoTrans.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         LDA >= max(1,M), if trans == Pasma[Conj]Trans.
 *         LDA >= max(1,N), if trans == PasmaNoTrans.
 *
 * @param[out] B
 *         Matrix of size LDB-by-M, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDB-by-N, if trans == PasmaNoTrans.
 *
 * @param[in] LDB
 *         The leading dimension of the array B.
 *         LDB >= max(1,N), if trans == Pasma[Conj]Trans.
 *         LDB >= max(1,M), if trans == PasmaNoTrans.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */
#if defined(CHAMELEON_HAVE_WEAK)
#pragma weak CORE_zlatro = PCORE_zlatro
#define CORE_zlatro PCORE_zlatro
#endif
int CORE_zlatro(cham_uplo_t uplo, cham_trans_t trans,
                int M, int N,
                const CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *B, int LDB)
{
    int i, j;

    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower) && (uplo != ChamUpperLower) ) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if ( !isValidTrans( trans ) ) {
        coreblas_error(2, "Illegal value of trans");
        return -2;
    }
    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ( (LDA < chameleon_max(1,M)) && (M > 0) ) {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < chameleon_max(1,N)) && (N > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    if (trans == ChamNoTrans) {
        CORE_zlacpy(uplo, M, N, A, LDA, B, LDB);
    }
    else {
#if defined(PRECISION_z) || defined(PRECISION_c)
        if (trans == ChamConjTrans)
        {
            if(uplo == ChamUpper) {
                for(j=0; j<N; j++)
                    for(i=0; i<chameleon_min(j+1,M); i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
            else if(uplo == ChamLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
        }
        else
#endif
        {
            if(uplo==ChamUpper) {
                for(j=0;j<N;j++)
                    for(i=0;i<chameleon_min(j+1,M);i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else if(uplo==ChamLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
        }
    }

    return CHAMELEON_SUCCESS;
}
