/**
 *
 * @file core_zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zgetrf with partial pivoting CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Extracts the rows from the tile A that have been selected as
 * pivot into the tile B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] m0
 *         The index of the first row of the tile A into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix A.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix B. This is the number of potential
 *         pivot that can be extracted from A.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to extract the
 *          pivot rows if some are selected in the range m0..m0+m.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n with 0s or already collected
 *          rows.
 *          On exit, B is filled with the selected rows from A, such that for
 *          each row i, B[i,:] = A[perm[i]-m0-1,:].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,k).
 *
 * @param[in] perm
 *          The permutation array of dimension k.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zlaswp_get( cham_side_t side, int m0, int m, int n, int k,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex64_t       *B, int ldb,
                 const int *perm )
{
    int i, idx_m, nb_elt;
    int A_ld, A_inc;
    int B_ld, B_inc;

    /* Check input arguments */
    if (m0 < 0) {
        coreblas_error(1, "Illegal value of m0");
        return -1;
    }
    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        coreblas_error(4, "Illegal value of k");
        return -4;
    }
    if (side == ChamLeft ) {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,k)) && (k > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }
    else {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( side == ChamLeft ) {
        A_inc  = 1;
        A_ld   = lda;
        B_inc  = 1;
        B_ld   = ldb;
        idx_m  = m;
        nb_elt = n;
    }
    else {
        A_inc  = lda;
        A_ld   = 1;
        B_inc  = ldb;
        B_ld   = 1;
        idx_m  = n;
        nb_elt = m;
    }

    for( i=0; i<k; i++ )
    {
        int idx = perm[i] - m0;

        if ( ( idx >= 0 ) && (idx < idx_m ) )
        {
            cblas_zcopy( nb_elt, A + idx * A_inc, A_ld,
                                 B + i   * B_inc, B_ld );
        }
    }

    return CHAMELEON_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Extracts the rows from the tile B that have been selected
 * as pivot into the tile A using the converted perm into set of couples.
 * See CORE_perm_to_idx()
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] m0
 *         The index of the first row/column of the tile A into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix A.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix B. This is the number of potential
 *         pivot that can be extracted from A.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to extract the
 *          pivot rows/columns if some are selected in the range m0..m0+m.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n with 0s or already collected
 *          rows.
 *          On exit, B is filled with the selected rows from A, such that for
 *          each row i, B[i] = A[perm[i]-m0-1].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,k).
 *
 * @param[in] perm_m
 *          The index of the tile in the perm_idx array.
 *
 * @param[in] perm_mt
 *          The number of tiles stored in perm_idx to compute the shift to the
 *          perm_data array.
 *
 * @param[in] perm_idx
 *          The permutation idx array of dimension (perm_mt+1) + 2 * m as defined in
 *          CORE_perm_to_idx().
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zlaswp_get_idx( cham_side_t side, int m0, int m, int n, int k,
                     const CHAMELEON_Complex64_t *A, int lda,
                     CHAMELEON_Complex64_t       *B, int ldb,
                     int perm_m, int perm_mt,
                     const int *perm_idx )
{
    int        i, nb_elt;
    int        A_ld, A_inc;
    int        B_ld, B_inc;
    int        first, last;
    const int *perm_data = perm_idx + perm_mt + 1;

    /* Check input arguments */
    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        coreblas_error(4, "Illegal value of k");
        return -4;
    }
    if (side == ChamLeft ) {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,k)) && (k > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }
    else {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( side == ChamLeft ) {
        A_inc  = 1;
        A_ld   = lda;
        B_inc  = 1;
        B_ld   = ldb;
        nb_elt = n;
    }
    else {
        A_inc  = lda;
        A_ld   = 1;
        B_inc  = ldb;
        B_ld   = 1;
        nb_elt = m;
    }

    first = perm_idx[ perm_m ];
    last  = perm_idx[ perm_m + 1 ];
    perm_data += first * 2;

    for( i = first; i < last; i++ ) {
        int ib = perm_data[ 0 ];
        int ia = perm_data[ 1 ] - m0;

        cblas_zcopy( nb_elt, A + ia * A_inc, A_ld,
                             B + ib * B_inc, B_ld );

        perm_data += 2;
    }

    return CHAMELEON_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Copies the rows from the tile A into the tile B when they are
 * the destination of the pivoted rows.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         Specifies the side of the permutation
 *         - ChamLeft  : rows permutation
 *         - ChamRight : columns permutation
 *
 * @param[in] m0
 *         The index of the first row of the tile B into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix B.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix A. This is the number of potential
 *         pivot that can be inserted into B.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to read the
 *          pivoted rows.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n that may require some pivoted rows.
 *          On exit, B is updated with the pivoted rows it needs to receive, such that for
 *          each row i, A[i,:] = B[invp[i]-m0-1,:].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[in] invp
 *          The inverse permutation array of dimension k.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zlaswp_set( cham_side_t side,
                 int m0, int m, int n, int k,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex64_t *B, int ldb,
                 const int *invp )
{
    int i, idx_m, nb_elt;
    int A_ld, A_inc;
    int B_ld, B_inc;

    /* Check input arguments */
    if (m0 < 0) {
        coreblas_error(1, "Illegal value of m0");
        return -1;
    }
    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        coreblas_error(4, "Illegal value of k");
        return -4;
    }

    if ( side == ChamLeft ) {
        if ((lda < chameleon_max(1,k)) && (k > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }
    else {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( side == ChamLeft ) {
        A_inc  = 1;
        A_ld   = lda;
        B_inc  = 1;
        B_ld   = ldb;
        idx_m  = m;
        nb_elt = n;
    }
    else {
        A_inc  = lda;
        A_ld   = 1;
        B_inc  = ldb;
        B_ld   = 1;
        idx_m  = n;
        nb_elt = m;
    }

    for( i=0; i<k; i++ )
    {
        int idx = invp[i] - m0;

        if ( ( idx >= 0 ) && (idx < idx_m ) ) {
            cblas_zcopy( nb_elt, A + i   * A_inc, A_ld,
                                 B + idx * B_inc, B_ld );
        }
    }

    return CHAMELEON_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * @brief Copies the rows from the tile A into the tile B when they are the
 * destination of the pivoted rows based on the invp_idx array as defined in
 * CORE_perm_to_idx().
 *
 *******************************************************************************
 *
 * @param[in] side
 *         Specifies the side of the permutation
 *         - ChamLeft  : rows permutation
 *         - ChamRight : columns permutation
 *
 * @param[in] m0
 *         The index of the first row of the tile B into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix B.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix A. This is the number of potential
 *         pivot that can be inserted into B.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to read the
 *          pivoted rows.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n that may require some pivoted rows.
 *          On exit, B is updated with the pivoted rows it needs to receive, such that for
 *          each row i, A[i] = B[invp[i]-m0-1].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[in] invp_m
 *          The index of the tile in the invp_idx array.
 *
 * @param[in] invp_mt
 *          The number of tiles stored in invp_idx to compute the shift to the
 *          invp_data array.
 *
 * @param[in] invp_idx
 *          The inverse permutation idx array of dimension (invp_mt+1) + 2 * m as
 *          defined in CORE_perm_to_idx().
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zlaswp_set_idx( cham_side_t side, int m0, int m, int n, int k,
                     const CHAMELEON_Complex64_t *A, int lda,
                     CHAMELEON_Complex64_t       *B, int ldb,
                     int invp_m, int invp_mt,
                     const int *invp_idx )
{
    int        i, nb_elt;
    int        A_ld, A_inc;
    int        B_ld, B_inc;
    int        first, last;
    const int *invp_data = invp_idx + invp_mt + 1;

    /* Check input arguments */
    if (m0 < 0) {
        coreblas_error(1, "Illegal value of m0");
        return -1;
    }
    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        coreblas_error(4, "Illegal value of k");
        return -4;
    }
    if (side == ChamLeft ) {
        if ((lda < chameleon_max(1,k)) && (k > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }
    else {
        if ((lda < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(6, "Illegal value of lda");
            return -6;
        }
        if ((ldb < chameleon_max(1,m)) && (m > 0)) {
            coreblas_error(8, "Illegal value of ldb");
            return -8;
        }
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    if ( side == ChamLeft ) {
        A_inc  = 1;
        A_ld   = lda;
        B_inc  = 1;
        B_ld   = ldb;
        nb_elt = n;
    }
    else {
        A_inc  = lda;
        A_ld   = 1;
        B_inc  = ldb;
        B_ld   = 1;
        nb_elt = m;
    }

    first = invp_idx[ invp_m ];
    last  = invp_idx[ invp_m + 1 ];
    invp_data += first * 2;

    for( i = first; i < last; i++ ) {
        int ia = invp_data[ 0 ];
        int ib = invp_data[ 1 ] - m0;

        cblas_zcopy( nb_elt, A + ia * A_inc, A_ld,
                             B + ib * B_inc, B_ld );

        invp_data += 2;
    }

    return CHAMELEON_SUCCESS;
}
