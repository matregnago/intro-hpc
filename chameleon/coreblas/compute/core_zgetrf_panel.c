/**
 *
 * @file core_zgetrf_panel.c
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
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

static const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;
static const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t)1.0;

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * CORE_zgetrf_panel_diag computes the LU factorization of a single column on
 * the diagonal block. Works in combination with CORE_zgetrf_panel_offdiag() to
 * perform the full LU factorization with partial pivoting of a single column in
 * a panel.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in] h
 *         The index of the column to factorize in the matrix A.
 *
 * @param[in] m0
 *         The number of rows above the diagonale tile A in the global matrix to
 *         be factorized.
 *
 * @param[in,out] A
 *          On entry, the matrix A where column h-1 needs to be factorized, and
 *          pivot for column h needs to be selected.
 *          The lower pentagonal part m-by-(h-1) are the L factors of the LU
 *          decomposition of A for the first (h-1) columns/rows.
 *          The upper pentagonal part (h-1)-by-n are the U factors of the LU
 *          decomposition of A for the first (h-1) columns/rows.
 *          On exit, A is factorized up to the (h-1)^th column/row and the
 *          remaining part is updated. The pivot of column h is selected.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[in,out] IPIV
 *          On entry, the pivot array of size min(m,n) with the first h-2
 *          columns initialized.
 *          On exit, IPIV[h-1] is updated with the selected pivot for the
 *          previous column.
 * *
 * @param[in,out] nextpiv
 *          On entry, the allocated and initialized CHAM_piv_t structure to
 *          store the information related to pivot at stage h.
 *          On exit, the diagrow is initialized with the diagonal h row from A,
 *          and the pivot row is updated if a better pivot is found for step h.
 *
 * @param[in] prevpiv
 *          Holds the CHAM_piv_t structure of the previous (h-1) step with the
 *          (h-1) diagonal row, and the pivot selected row for (h-1) step.
 *          Not referenced if h is equal to 0.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zgetrf_panel_diag( int m, int n, int h, int m0, int ib,
                        CHAMELEON_Complex64_t *A, int lda,
                        CHAMELEON_Complex64_t *U, int ldu,
                        int *IPIV, CHAM_pivot_t *nextpiv, const CHAM_pivot_t *prevpiv )
{
    CHAMELEON_Complex64_t *subA = A + lda * h + h;
    CHAMELEON_Complex64_t  alpha;

    /* Compute localn, the number of columns to update inside the block */
    int localnt = chameleon_ceil( chameleon_min(m, n), ib );
    int localn  = (( h / ib ) == (localnt-1)) ? n - ( h / ib ) * ib : ib;

    if ( h > 0 ) {
        CHAMELEON_Complex64_t *L, *pivrow;

        int pivot_blkm0  = prevpiv->blkm0;
        int pivot_blkidx = prevpiv->blkidx;
        int pivot        = pivot_blkm0 + pivot_blkidx;

        pivrow = (CHAMELEON_Complex64_t *)(prevpiv->pivrow);

        /*
         * Copy the previously selected row to the right place if it's not the original
         * diagonal row or if does not belong to the diagonal block
         */
        if(   (pivot_blkm0 != m0) ||
            ( (pivot_blkm0 == m0) && (pivot_blkidx != (h-1)) ) )
        {
            cblas_zcopy( n, pivrow, 1,
                         A + h-1, lda );
        }

        /* Store the pivot */
        IPIV[h-1] = pivot + 1;

        /*
         * Copy the diagonal row in place of the selected one (previous pivot) to the right place if:
         *    - we own the previous pivot
         *    - pivot index was not last step index
         */
        if( (pivot_blkm0  == m0 ) &&
            (pivot_blkidx != h-1) )
        {
            cblas_zcopy( n, prevpiv->diagrow, 1,
                            A + pivot_blkidx, lda );
        }

        /* Take the pointer to the element below the current diagonal element */
        L = A + lda * (h-1) + h;

        if ( h < m ) {
            /* Scaling off-diagonal terms */
            alpha = (CHAMELEON_Complex64_t)1. / pivrow[h-1];
            cblas_zscal( m-h, CBLAS_SADDR( alpha ), L, 1 );
        }

        if ( h < chameleon_min( m, n ) ) {
            /* Applying the update */
            if ( (h % ib) != 0 ) {
                /* We apply a BLAS-2 operation with a single outer-product */
                cblas_zgeru( CblasColMajor, m-h, localn-(h%ib),
                             CBLAS_SADDR(mzone),
                             L,          1,
                             pivrow + h, 1,
                             subA,       lda );
            }
            else {
                /* We apply a BLAS-3 operation with a matrix-matrix product */
                const CHAMELEON_Complex64_t *L    = A + (h-ib) * lda + h;
                const CHAMELEON_Complex64_t *subU = U +  h     * ldu;

                LAPACKE_zlacpy_work(
                    LAPACK_COL_MAJOR, 'A', ib, n - h,
                    subU, ldu, A + h * lda + (h-ib), lda );

                cblas_zgemm( CblasColMajor,
                             (CBLAS_TRANSPOSE)CblasNoTrans, (CBLAS_TRANSPOSE)CblasNoTrans,
                             m - h, n - h, ib,
                             CBLAS_SADDR(mzone), L,    lda,
                                                 subU, ldu,
                             CBLAS_SADDR( zone), subA, lda );
            }
        }
        else {
          return 0;
        }
    }

    /* Looking for the local max index */
    {
        CHAMELEON_Complex64_t *next_pivrow = (CHAMELEON_Complex64_t *)(nextpiv->pivrow);
        int                    index;

        index = cblas_izamax( m-h, subA, 1 );
        alpha = subA[index];

        if ( cabs(alpha) > cabs(next_pivrow[h]) ) {
            /* We save this index as a potential pivot for the next column of the panel  */
            nextpiv->blkm0  = m0;                   /* block index */
            nextpiv->blkidx = index + h;            /* index in local numbering */

            /* Save the row holding the selected pivot */
            cblas_zcopy( n, A + index + h, lda,
                            next_pivrow,   1 );
        }
    }

    /* Store current diagonal row (in full) into pivot structure */
    cblas_zcopy( n, A + h,            lda,
                    nextpiv->diagrow, 1 );

    return 0;
}

/**
 ******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 * CORE_zgetrf_panel_offdiag computes the LU factorization of a single column on
 * an off-diagonal block. Works in combination with CORE_zgetrf_panel_diag() to
 * perform the full LU factorization with partial pivoting of a single column in
 * a panel.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in] h
 *         The index of the column to factorize in the matrix A.
 *
 * @param[in,out] A
 *          On entry, the matrix A where column h-1 needs to be factorized, and
 *          pivot for column h needs to be selected.
 *          The lower pentagonal part m-by-(h-1) are the L factors of the LU
 *          decomposition of A for the first (h-1) columns/rows.
 *          The upper pentagonal part (h-1)-by-n are the U factors of the LU
 *          decomposition of A for the first (h-1) columns/rows.
 *          On exit, A is factorized up to the (h-1)^th column/row and the
 *          remaining part is updated. The pivot of column h is selected.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[in,out] nextpiv
 *          On entry, the allocated and initialized CHAM_piv_t structure to
 *          store the information related to pivot at stage h.
 *          On exit, the diagrow is initialized with the diagonal h row from A,
 *          and the pivot row is updated if a better pivot is found for step h.
 *
 * @param[in] prevpiv
 *          Holds the CHAM_piv_t structure of the previous (h-1) step with the
 *          (h-1) diagonal row, and the pivot selected row for (h-1) step.
 *          Not referenced if h is equal to 0.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zgetrf_panel_offdiag( int m, int n, int h, int m0, int ib,
                           CHAMELEON_Complex64_t *A, int lda,
                           CHAMELEON_Complex64_t *U, int ldu,
                           CHAM_pivot_t *nextpiv, const CHAM_pivot_t *prevpiv )
{
    CHAMELEON_Complex64_t *subA = A + lda * h;
    CHAMELEON_Complex64_t  alpha;

    /* Compute the number of columns to update inside the block */
    int localnt = chameleon_ceil( n, ib );
    int localn  = (( h / ib ) == (localnt-1)) ? n - ( h / ib ) * ib : ib;

    if ( h != 0 ) {
        CHAMELEON_Complex64_t *L, *pivrow;

        int pivot_blkm0  = prevpiv->blkm0;
        int pivot_blkidx = prevpiv->blkidx;

        pivrow = (CHAMELEON_Complex64_t *)(prevpiv->pivrow);

        /*
         * Copy the diagonal row in place of the selected one (previous pivot) to the right place if:
         *    - the previous pivot block was this one
         *    - it was not the previous row ??
         */
        if( pivot_blkm0 == m0 )
        {
            cblas_zcopy( n, prevpiv->diagrow, 1,
                            A + pivot_blkidx, lda );
        }

        /* Take the pointer to the element below the current diagonal element */
        L = A + lda * (h-1);

        /* Scaling off-diagonal terms */
        alpha = (CHAMELEON_Complex64_t)1. / pivrow[h-1];
        cblas_zscal( m, CBLAS_SADDR( alpha ), L, 1 );

        /*
         * h is compared only to n, because if we are on the last column of a
         * tile, m might be much smaller than n, and still we need to apply
         * the geru call. If this is the diagonal tile, we will just look for
         * the next maximum for nothing.
         */
        if ( h < n ) {
            /* Applying the update */
            if ( (h % ib) != 0 ) {
                /* We apply a BLAS-2 operation with a single outer-product */
                cblas_zgeru( CblasColMajor, m, localn-(h%ib),
                             CBLAS_SADDR(mzone),
                             L,          1,
                             pivrow + h, 1,
                             subA,       lda );
            }
            else {
                /* We apply a BLAS-3 operation with a matrix-matrix product */
                const CHAMELEON_Complex64_t *L    = A + (h-ib) * lda;
                const CHAMELEON_Complex64_t *subU = U +  h     * ldu;

                cblas_zgemm( CblasColMajor,
                             (CBLAS_TRANSPOSE)CblasNoTrans, (CBLAS_TRANSPOSE)CblasNoTrans,
                             m, n - h, ib,
                             CBLAS_SADDR(mzone), L,    lda,
                                                 subU, ldu,
                             CBLAS_SADDR( zone), subA, lda );
            }
        }
        else {
            return 0;
        }
    }

    /* Looking for the local max index */
    {
        CHAMELEON_Complex64_t *next_pivrow = (CHAMELEON_Complex64_t *) nextpiv->pivrow;
        int                    index;

        index = cblas_izamax( m, subA, 1 );
        alpha = subA[index];

        if ( cabs(alpha) > cabs(next_pivrow[h]) )
        {
            /* We save this index as a potential pivot for the next column of the panel */
            nextpiv->blkm0  = m0;           /* block index */
            nextpiv->blkidx = index;        /* index in local numbering */

            /* Save the row holding the selected pivot */
            cblas_zcopy( n, A + index,   lda,
                            next_pivrow, 1 );
        }
    }
    return 0;
}
