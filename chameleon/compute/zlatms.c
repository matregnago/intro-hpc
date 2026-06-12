/**
 *
 * @file zlatms.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatms wrappers
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Generates random matrices with specified singular values
 *    (or hermitian with specified eigenvalues)
 *    for testing programs (See LAPACK).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0. When M = 0,
 *          the returned value is set to zero.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] idist
 *          On entry, idist specifies the type of distribution to be used
 *          to generate the random eigen-/singular values.
 *          ChamDistUniform   => UNIFORM( 0, 1 )  ( 'U' for uniform )
 *          ChamDistSymmetric => UNIFORM( -1, 1 ) ( 'S' for symmetric )
 *          ChamDistNormal    => NORMAL( 0, 1 )   ( 'N' for normal )
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 * @param[in] sym
 *          - ChamHermGeev, the generated matrix is hermitian, with
 *            eigenvalues specified by D, cond, mode, and dmax; they
 *            may be positive, negative, or zero.
 *          - ChamHermPoev, the generated matrix is hermitian, with
 *            eigenvalues (= singular values) specified by D, cond,
 *            mode, and dmax; they will not be negative.
 *          - ChamNonsymPosv, the generated matrix is nonsymmetric, with
 *            singular values specified by D, cond, mode, and dmax;
 *            they will not be negative.
 *          - ChamSymPosv, the generated matrix is (complex) symmetric,
 *            with singular values specified by D, cond, mode, and
 *            dmax; they will not be negative.
 *
 * @param[in,out] D
 *          Array of dimension ( MIN( M, N ) )
 *          This array is used to specify the singular values or
 *          eigenvalues of A (see SYM, above.)  If mode=0, then D is
 *          assumed to contain the singular/eigenvalues, otherwise
 *          they will be computed according to mode, cond, and dmax,
 *          and placed in D.
 *          Modified if mode is nonzero.
 *
 * @param[in] mode
 *          On entry this describes how the singular/eigenvalues are to
 *          be specified:
 *          mode = 0 means use D as input
 *          mode = 1 sets D(1)=1 and D(2:N)=1.0/cond
 *          mode = 2 sets D(1:N-1)=1 and D(N)=1.0/cond
 *          mode = 3 sets D(I)=cond**(-(I-1)/(N-1))
 *          mode = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/cond)
 *          mode = 5 sets D to random numbers in the range
 *                   ( 1/cond , 1 ) such that their logarithms
 *                   are uniformly distributed.
 *          mode = 6 set D to random numbers from same distribution
 *                   as the rest of the matrix.
 *          mode < 0 has the same meaning as ABS(mode), except that
 *             the order of the elements of D is reversed.
 *          Thus if mode is positive, D has entries ranging from
 *             1 to 1/cond, if negative, from 1/cond to 1,
 *          If sym == ChamHermPoev, and mode is neither 0, 6, nor -6, then
 *             the elements of D will also be multiplied by a random
 *             sign (i.e., +1 or -1.)
 *
 * @param[in] cond
 *          On entry, this is used as described under mode above.
 *          If used, it must be >= 1.
 *
 * @param[in] dmax
 *          If mode is neither -6, 0 nor 6, the contents of D, as
 *          computed according to mode and cond, will be scaled by
 *          dmax / max(abs(D(i))); thus, the maximum absolute eigen- or
 *          singular value (which is to say the norm) will be abs(dmax).
 *          Note that dmax need not be positive: if dmax is negative
 *          (or zero), D will be scaled by a negative number (or zero).
 *          Not modified.
 *
 * @param[out] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval the two-norm estimate.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatms_Tile
 * @sa CHAMELEON_zlatms_Tile_Async
 * @sa CHAMELEON_clange
 * @sa CHAMELEON_dlange
 * @sa CHAMELEON_slange
 *
 */
int CHAMELEON_zlatms( int M, int N, cham_dist_t idist,
                      unsigned long long int seed, cham_sym_t sym,
                      double *D, int mode, double cond, double dmax,
                      CHAMELEON_Complex64_t *A, int LDA )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlatms", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescOutput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zlatms_Tile_Async( idist, seed, sym, D, mode, cond, dmax,
                                 &descAt, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescOutput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Tile equivalent of CHAMELEON_zlatms().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On entry, the input matrix A.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatms
 * @sa CHAMELEON_zlatms_Tile_Async
 * @sa CHAMELEON_clange_Tile
 * @sa CHAMELEON_dlange_Tile
 * @sa CHAMELEON_slange_Tile
 *
 */
int CHAMELEON_zlatms_Tile( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym,
                           double *D, int mode, double cond, double dmax, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlatms_Tile_Async( idist, seed, sym, D, mode, cond, dmax, A, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 * @brief Non-blocking equivalent of CHAMELEON_zlatms_Tile().
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatms
 * @sa CHAMELEON_zlatms_Tile
 * @sa CHAMELEON_clange_Tile_Async
 * @sa CHAMELEON_dlange_Tile_Async
 * @sa CHAMELEON_slange_Tile_Async
 *
 */
int CHAMELEON_zlatms_Tile_Async( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym,
                                 double *D, int mode, double cond, double dmax, CHAM_desc_t *A,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms_Tile", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlatms_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_zlatms_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min(A->m, A->n) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzlatms( idist, seed, sym, D, mode, cond, dmax, A, sequence, request );

    return CHAMELEON_SUCCESS;
}
