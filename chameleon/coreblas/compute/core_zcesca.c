/**
 *
 * @file core_zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zcesca CPU kernel
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zcesca computes the Centered-Scaled matrix factors of A inplace
 *
 *  Considering a matrix A of size m x n, \f[A = (a_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f]
 *  Lets
 *  \f[g_i = \frac{1}{n} \sum_j a_{ij} \\
 *     g_j = \frac{1}{m} \sum_i a_{ij} \\
 *     g   = \frac{1}{mn} \sum_{i,j} a_{ij}\f]
 *  A centered rowwise gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_i \f]
 *  A centered columnwise gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_j \f]
 *  A bicentered gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 *  \f[ \bar{a}_{i,j} = a_{i,j} - g_i - g_j + g \f]
 * Lets
 * \f[d_i = || a_{i*} || = \sqrt{ \sum_j a_{ij}^2} \\
 *    d_j = || a_{*j} || = \sqrt{ \sum_i a_{ij}^2} \f]
 * A scaled rowwise gives \f[A' = (a_{i,j}')_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ a_{i*}' = \frac{a_{i*}}{d_i} \f]
 * A scaled columnwise gives \f[A' = (a_{i,j}')_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ a_{*j}' = \frac{a_{*j}}{d_j} \f]
 *
 * This function can also be used to compute a pretreatment of a Correspondence analysis (COA).
 * To use it set center = 1, scale = 1, axis = ChamEltwise.
 * In this case lets
 * \f[r_i = \sum_j a_{ij} \\
 *    c_j = \sum_i a_{ij} \\
 *    sg  = \sum_{i,j} a_{ij}\f]
 * A transformed gives \f[\bar{A} = (\bar{a}_{i,j})_{1 \leq i \leq m, 1 \leq j \leq n}\f] such that
 * \f[ \bar{a}_{i,j} = \frac{a_{i,j}-r_i*c_j/sg}{ \sqrt{r_i*c_j} } \f]
 *
 *******************************************************************************
 *
 * @param[in] center
 *          1 if A must be centered, else 0.
 *
 * @param[in] scale
 *          1 if A must be scaled, else 0.
 *
 * @param[in] axis
 *          Specifies the axis over which to center and or scale.
 *            = ChamColumnwise: centered column-wise
 *            = ChamRowwise: centered row-wise
 *            = ChamEltwise:
 *              bi-centered if center=1 and scale=0
 *              COA if center=1 and scale=1
 *
 * @param[in] M
 *          The number of rows of the overall matrix.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the overall matrix.  N >= 0.
 *
 * @param[in] Mt
 *          The number of rows of the tile A.  Mt >= 0.
 *
 * @param[in] Nt
 *         The number of columns of the tile A.  Nt >= 0.
 *
 * @param[in] Gi
 *         The 1-by-Mt tile containing the sum of values by rows.
 *
 * @param[in] LDGI
 *         The leading dimension of the array Gi.
 *
 * @param[in] Gj
 *         The 1-by-Nt tile containing the sum of values by columns.
 *
 * @param[in] LDGJ
 *         The leading dimension of the array Gj.
 *
 * @param[in] G
 *         The sum of values of all the matrix.
 *
  * @param[in] Di
 *         The 1-by-Mt tile containing the sum of squares by rows.
 *
 * @param[in] LDDI
 *         The leading dimension of the array Di.
 *
 * @param[in] Dj
 *         The 1-by-Nt tile containing the sum of squares by columns.
 *
 * @param[in] LDDJ
 *         The leading dimension of the array Dj.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the centered-scaled matrix coefficients.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 */

int CORE_zcesca( int center, int scale,
                 cham_store_t axis,
                 int M, int N,
                 int Mt, int Nt,
                 const CHAMELEON_Complex64_t *Gi, int LDGI,
                 const CHAMELEON_Complex64_t *Gj, int LDGJ,
                 const CHAMELEON_Complex64_t *G,
                 const double *Di, int LDDI,
                 const double *Dj, int LDDJ,
                 CHAMELEON_Complex64_t *A, int LDA )
{
    int i, j;

    /* Check input arguments */
    if ( (center != 0) && (center != 1) ) {
        coreblas_error(1, "Illegal value of center");
        return -1;
    }
    if ( (scale != 0) && (scale != 1) ) {
        coreblas_error(2, "Illegal value of scale");
        return -2;
    }
    if ( (axis != ChamColumnwise) && (axis != ChamRowwise) && (axis != ChamEltwise) ) {
        coreblas_error(3, "Illegal value of axis");
        return -3;
    }
    if (M < 0) {
        coreblas_error(4, "Illegal value of M");
        return -4;
    }
    if (N < 0) {
        coreblas_error(5, "Illegal value of N");
        return -5;
    }
    if (Mt < 0) {
        coreblas_error(6, "Illegal value of Mt");
        return -6;
    }
    if (Nt < 0) {
        coreblas_error(7, "Illegal value of Nt");
        return -7;
    }
    if (LDGI < 0) {
        coreblas_error(9, "Illegal value of LDGI");
        return -9;
    }
    if (LDGJ < 0) {
        coreblas_error(11, "Illegal value of LDGJ");
        return -11;
    }
    if (LDDI < 0) {
        coreblas_error(14, "Illegal value of LDDI");
        return -14;
    }
    if (LDDJ < 0) {
        coreblas_error(16, "Illegal value of LDDJ");
        return -16;
    }
    if ((LDA < chameleon_max(1,Mt)) && (Mt > 0)) {
        coreblas_error(18, "Illegal value of LDA");
        return -18;
    }

    /* Quick return */
    if ((center == 0) && (scale == 0)) {
        return CHAMELEON_SUCCESS;
    }
    if ((Mt == 0) || (Nt == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( !( (center == 1) && (scale == 1) && (axis == ChamEltwise) ) ) {
        /* PCA case i.e. centered-scaled or bi-centering */
        CHAMELEON_Complex64_t gi = (CHAMELEON_Complex64_t)0.;
        CHAMELEON_Complex64_t gj = (CHAMELEON_Complex64_t)0.;
        CHAMELEON_Complex64_t g  = (CHAMELEON_Complex64_t)0.;
        double di = 0.;
        double dj = 0.;

        if ( (center == 1) && (axis == ChamEltwise) ) {
            /* overall mean of values */
            g =  G[0] / ( (double)M * (double)N );
        }

        for(j = 0; j < Nt; j++) {
            if ( (center == 1) && ( (axis == ChamColumnwise) || (axis == ChamEltwise) ) ) {
                /* mean of values of the column */
                gj = Gj[j*LDGJ] / ((double)M);
            }
            if ( (scale == 1) && (axis == ChamColumnwise) ) {
                /* norm 2 of the column */
                dj = Dj[j*LDDJ];
            }
            for(i = 0; i < Mt; i++) {
                if ( center == 1 ) {
                    if ( (axis == ChamRowwise) || (axis == ChamEltwise) ) {
                        /* mean of values of the row */
                        gi = Gi[i] / ((double)N);
                        /* compute centered matrix factor */
                        A[j*LDA+i] -= gi;
                    }
                    if ( (axis == ChamColumnwise) || (axis == ChamEltwise) ) {
                        /* compute centered matrix factor */
                        A[j*LDA+i] -= gj;
                    }
                    if ( axis == ChamEltwise ) {
                        /* compute centered matrix factor */
                        A[j*LDA+i] += g;
                    }
                }
                if ( scale == 1 ) {
                    if ( axis == ChamColumnwise ) {
                        /* compute scaled matrix factor */
                        A[j*LDA+i] /= dj;
                    }
                    if ( axis == ChamRowwise ) {
                        /* norm 2 of the row */
                        di = Di[i];
                        /* compute scaled matrix factor */
                        A[j*LDA+i] /= di;
                    }
                }
            }
        }

    } else {
        /* COA case */
        CHAMELEON_Complex64_t rc, sqrc;

        /* update the matrix */
        for(j = 0; j < Nt; j++) {
            for(i = 0; i < Mt; i++) {
                rc = Gi[i] * Gj[j];
                sqrc = sqrt(rc);
                A[j*LDA + i] -= Gi[i] * Gj[j] / G[0];
                A[j*LDA + i] /= sqrc;
            }
        }

    }

    return CHAMELEON_SUCCESS;
}
