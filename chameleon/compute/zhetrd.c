/**
 *
 * @file zhetrd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhetrd wrappers
 *
 * @version 1.4.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#endif

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zhetrd - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary
 *  Q2).  Let Q = Q1 * Q2 be the global unitary transformation; Q^H *
 *  A * Q = S.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = ChamNoVec: computes tridiagonal only;
 *          = ChamVec: computes tridiagonal and generate the orthogonal matrix Q.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = ChamUpper:: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = ChamLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = ChamLower) or the
 *          upper triangle (if uplo = ChamUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = ChamUpper, E(i) = A(i+1,i) if uplo = ChamLower.
 *
 * @param[out] descT
 *          On entry, descriptor as return by CHAMELEON_Alloc_Workspace_zhetrd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = ChamVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = ChamNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhetrd_Tile
 * @sa CHAMELEON_zhetrd_Tile_Async
 * @sa CHAMELEON_chetrd
 * @sa CHAMELEON_dsytrd
 * @sa CHAMELEON_ssytrd
 *
 */
int CHAMELEON_zhetrd( cham_job_t jobz, cham_uplo_t uplo, int N,
                      CHAMELEON_Complex64_t *A, int LDA,
                      double *D,
                      double *E,
                      CHAM_desc_t *descT,
                      CHAMELEON_Complex64_t *Q, int LDQ )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_zhetrd", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != ChamNoVec && jobz != ChamVec) {
        chameleon_error("CHAMELEON_zhetrd", "illegal value of jobz");
        return -1;
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zhetrd", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zhetrd", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAMELEON_zhetrd", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (N == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZHETRD, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhetrd", "chameleon_tune() failed");
        return status;
    }
    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, uplo,
                         A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zhetrd_Tile_Async( jobz, uplo, &descAt, D, E, descT, Q, LDQ, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zhetrd_Tile - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary Q2).
 *  Let Q = Q1 * Q2 be the global unitary transformation;
 *  Q^H * A * Q = S.
 *  Tile equivalent of CHAMELEON_zhetrd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = ChamNoVec: computes tridiagonal only;
 *          = ChamVec: computes tridiagonal and generate the orthogonal matrix Q.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.  If uplo
 *          = ChamUpper, the leading N-by-N upper triangular part of
 *          A contains the upper triangular part of the matrix A, and
 *          the strictly lower triangular part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = ChamVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = ChamNoVec, then on exit the lower triangle (if
 *          uplo = ChamLower) or the upper triangle (if uplo =
 *          ChamUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = ChamUpper,
 *          E(i) = A(i+1,i) if uplo = ChamLower.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = ChamVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = ChamNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zhetrd
 * @sa CHAMELEON_zhetrd_Tile_Async
 * @sa CHAMELEON_chetrd_Tile
 * @sa CHAMELEON_dsytrd_Tile
 * @sa CHAMELEON_ssytrd_Tile
 * @sa CHAMELEON_zhetrd_Tile
 *
 */
int CHAMELEON_zhetrd_Tile( cham_job_t jobz, cham_uplo_t uplo,
                           CHAM_desc_t *A, double *D, double *E,
                           CHAM_desc_t *T, CHAMELEON_Complex64_t *Q, int LDQ )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhetrd_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zhetrd_Tile_Async( jobz, uplo, A, D, E, T, Q, LDQ, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );

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
 *  CHAMELEON_zhetrd_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
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
 * @sa CHAMELEON_zhetrd
 * @sa CHAMELEON_zhetrd_Tile
 * @sa CHAMELEON_chetrd_Tile_Async
 * @sa CHAMELEON_dsytrd_Tile_Async
 * @sa CHAMELEON_ssytrd_Tile_Async
 *
 */
int CHAMELEON_zhetrd_Tile_Async( cham_job_t jobz,
                                 cham_uplo_t uplo,
                                 CHAM_desc_t *A,
                                 double *W,
                                 double *E,
                                 CHAM_desc_t *T,
                                 CHAMELEON_Complex64_t *Q, int LDQ,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t descA;
    CHAM_desc_t descAB;
    int N, NB, LDAB;
    CHAM_desc_t D, *Dptr = NULL;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zhetrd_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zhetrd_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zhetrd_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (chameleon_desc_check(T) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check input arguments */
    if (jobz != ChamNoVec && jobz != ChamVec) {
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "illegal value of jobz");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if ((uplo != ChamLower) && (uplo != ChamUpper)) {
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "matrix need to be square");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        chameleon_error("CHAMELEON_zhetrd_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    N  = descA.m;
    NB = descA.mb;
#if defined(CHAMELEON_COPY_DIAG)
    {
        chameleon_ddesc_alloc_diag( chamctxt, &D, "HETRD_D", A->mb, A->m, A->n,
                                    chameleon_desc_datadist_get_iparam(A, 0),
                                    chameleon_desc_datadist_get_iparam(A, 1) );
        Dptr = &D;
    }
#endif
    /* Reduction to band. On exit, T contains reflectors */
    chameleon_pzhetrd_he2hb( uplo, A, T, Dptr,
                             sequence, request );

    LDAB = NB+1;

    /* Allocate band structure */
    chameleon_zdesc_alloc( descAB, "HETRD_AB",
                           LDAB, NB, /* mb, nb */
                           LDAB, N,  /* m, n */
                            );

    /* Copy data into band structure */
    chameleon_pztile2band( uplo, A, &descAB,
                           sequence, request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Reduce band matrix to tridiagonal matrix */
#if !defined(CHAMELEON_SIMULATION)
    {
        int info = LAPACKE_zhbtrd( LAPACK_COL_MAJOR,
                                   chameleon_lapack_const(jobz),
                                   chameleon_lapack_const(uplo),
                                   N, NB,
                                   (CHAMELEON_Complex64_t *) descAB.mat, LDAB,
                                   W, E, Q, LDQ );
        if (info != 0) {
            chameleon_error("CHAMELEON_zhetrd_Tile_Async", "LAPACKE_zhbtrd failed");
        }
    }
#endif /* !defined(CHAMELEON_SIMULATION) */
    if (Dptr != NULL) {
        chameleon_desc_destroy( Dptr );
    }
    chameleon_desc_destroy( &descAB );
    (void)D;
    return CHAMELEON_SUCCESS;
}
