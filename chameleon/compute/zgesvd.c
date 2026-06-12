/**
 *
 * @file zgesvd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesvd wrappers
 *
 * @version 1.4.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#endif

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgesvd - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *
 *       A = U * SIGMA * transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns V^T, not V.
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = ChamAllVec = 'A'(lapack): all M columns of U are returned
 *                         in array U;
 *          = ChamNoVec  = 'N': no columns of U (no left singular vectors)
 *                         are computed.
 *          = ChamSVec   = 'S': the first min(m,n) columns of U (the left
 *                         singular vectors) are returned in the array U;
 *                         NOT SUPPORTED YET
 *          = ChamOVec   = 'O': the first min(m,n) columns of U (the left
 *                         singular vectors) are overwritten on the array A;
 *                         NOT SUPPORTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V^H.
 *          Intended usage:
 *          = ChamAllVec = 'A'(lapack): all N rows of V^H are returned
 *                         in the array VT;
 *          = ChamNoVec  = 'N': no rows of V^H (no right singular vectors)
 *                         are computed.
 *          = ChamSVec   = 'S': the first min(m,n) rows of V^H (the right
 *                         singular vectors) are returned in the array VT;
 *                         NOT SUPPORTED YET
 *          = ChamOVec   = 'O': the first min(m,n) rows of V^H (the right
 *                         singular vectors) are overwritten on the array A;
 *                         NOT SUPPORTED YET
 *
 *          Note: jobu and jobvt cannot both be ChamOVec.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if jobu == ChamOVec, A is overwritten with the first min(m,n)
 *                               columns of U (the left singular vectors,
 *                               stored columnwise);
 *          if jobvt == ChamOVec, A is overwritten with the first min(m,n)
 *                                rows of V^H (the right singular vectors,
 *                                stored rowwise);
 *          if jobu != ChamOVec and jobvt != ChamOVec, the content of A
 *                                                     is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] S
 *          The double precision singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by CHAMELEON_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          If jobu == ChamAllVec, U contains the M-by-M unitary matrix U;
 *          if jobu == ChamSVec, U contains the first min(m,n) columns of U
 *                               (the left singular vectors, stored columnwise);
 *          if jobu == ChamNoVec or ChamOVec, U is not referenced.

 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1;
 *          if jobu == ChamSVec or ChamAllVec, LDU >= M.
 *
 * @param[out] VT
 *          If jobvt == ChamAllVec, VT contains the N-by-N unitary matrix
 *                                  V^H;
 *          if jobvt == ChamSVec, VT contains the first min(m,n) rows of
 *                                V^H (the right singular vectors, stored rowwise);
 *          if jobvt == ChamNoVec or ChamOVec, VT is not referenced.

 *
 * @param[in] LDVT
 *          The leading dimension of the array VT.  LDVT >= 1;
 *          if jobvt == ChamAllVec, LDVT >= N;
 *          if jobvt == ChamSVec,   LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesvd_Tile
 * @sa CHAMELEON_zgesvd_Tile_Async
 * @sa CHAMELEON_cgesvd
 * @sa CHAMELEON_dgesvd
 * @sa CHAMELEON_sgesvd
 *
 */
int CHAMELEON_zgesvd( cham_job_t jobu, cham_job_t jobvt,
                      int M, int N,
                      CHAMELEON_Complex64_t *A, int LDA,
                      double *S,
                      CHAM_desc_t *descT,
                      CHAMELEON_Complex64_t *U, int LDU,
                      CHAMELEON_Complex64_t *VT, int LDVT )
{
    int NB;
    int status;
    CHAM_context_t  *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgesvd", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    assert( chamctxt->scheduler != RUNTIME_SCHED_PARSEC );

    /* Check input arguments */
    if ( (jobu != ChamNoVec) && (jobu != ChamAllVec) ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of jobu" );
        return -1;
    }
    if ( (jobvt != ChamNoVec) && (jobvt != ChamAllVec) ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of jobvt" );
        return -2;
    }
    if ( M < 0 ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of M") ;
        return -3;
    }
    if ( N < 0 ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of N" );
        return -4;
    }
    if ( LDA < chameleon_max(1, M) ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of LDA" );
        return -6;
    }
    if ( LDU < 1 ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of LDU" );
        return -9;
    }
    if ( LDVT < 1 ) {
        chameleon_error( "CHAMELEON_zgesvd", "illegal value of LDVT" );
        return -11;
    }
    /* Quick return */
    if ( chameleon_min(M, N ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = chameleon_tune( CHAMELEON_FUNC_ZGESVD, M, N, 0 );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zgesvd", "chameleon_tune() failed" );
        return status;
    }

    /* Set MT, NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescInout, ChamUpperLower,
                         A, NB, NB,  LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    CHAMELEON_zgesvd_Tile_Async( jobu, jobvt, &descAt, S, descT, U, LDU, VT, LDVT, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, ChamUpperLower, sequence, &request );

    CHAMELEON_Desc_Flush( descT, sequence );
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
 *  CHAMELEON_zgesvd_Tile - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Tile equivalent of CHAMELEON_zgesvd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = ChamAllVec = 'A'(lapack): all M columns of U are returned
 *                         in array U;
 *          = ChamNoVec  = 'N': no columns of U (no left singular vectors)
 *                         are computed.
 *          = ChamSVec   = 'S': the first min(m,n) columns of U (the left
 *                         singular vectors) are returned in the array U;
 *                         NOT SUPPORTED YET
 *          = ChamOVec   = 'O': the first min(m,n) columns of U (the left
 *                         singular vectors) are overwritten on the array A;
 *                         NOT SUPPORTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V^H.
 *          Intended usage:
 *          = ChamAllVec = 'A'(lapack): all N rows of V^H are returned
 *                         in the array VT;
 *          = ChamNoVec  = 'N': no rows of V^H (no right singular vectors)
 *                         are computed.
 *          = ChamSVec   = 'S': the first min(m,n) rows of V^H (the right
 *                         singular vectors) are returned in the array VT;
 *                         NOT SUPPORTED YET
 *          = ChamOVec   = 'O': the first min(m,n) rows of V^H (the right
 *                         singular vectors) are overwritten on the array A;
 *                         NOT SUPPORTED YET
 *
 *          Note: jobu and jobvt cannot both be ChamOVec.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if jobu == ChamOVec, A is overwritten with the first min(m,n)
 *                               columns of U (the left singular vectors,
 *                               stored columnwise);
 *          if jobvt == ChamOVec, A is overwritten with the first min(m,n)
 *                                rows of V^H (the right singular vectors,
 *                                stored rowwise);
 *          if jobu != ChamOVec and jobvt != ChamOVec, the content of A
 *                                                     is destroyed.
 *
 * @param[out] S
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] T
 *          On entry, descriptor as return by CHAMELEON_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          If jobu == ChamAllVec, U contains the M-by-M unitary matrix U;
 *          if jobu == ChamSVec, U contains the first min(m,n) columns of U
 *                               (the left singular vectors, stored columnwise);
 *          if jobu == ChamNoVec or ChamOVec, U is not referenced.

 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1;
 *          if jobu == ChamSVec or ChamAllVec, LDU >= M.
 *
 * @param[out] VT
 *          If jobvt == ChamAllVec, VT contains the N-by-N unitary matrix
 *                                  V^H;
 *          if jobvt == ChamSVec, VT contains the first min(m,n) rows of
 *                                V^H (the right singular vectors, stored rowwise);
 *          if jobvt == ChamNoVec or ChamOVec, VT is not referenced.

 *
 * @param[in] LDVT
 *          The leading dimension of the array VT.  LDVT >= 1;
 *          if jobvt == ChamAllVec, LDVT >= N;
 *          if jobvt == ChamSVec,   LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgesvd
 * @sa CHAMELEON_zgesvd_Tile_Async
 * @sa CHAMELEON_cgesvd_Tile
 * @sa CHAMELEON_dgesvd_Tile
 * @sa CHAMELEON_sgesvd_Tile
 *
 */
int CHAMELEON_zgesvd_Tile( cham_job_t jobu, cham_job_t jobvt,
                           CHAM_desc_t *A,
                           double *S,
                           CHAM_desc_t *T,
                           CHAMELEON_Complex64_t *U, int LDU,
                           CHAMELEON_Complex64_t *VT, int LDVT )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgesvd_Tile", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    assert( chamctxt->scheduler != RUNTIME_SCHED_PARSEC );
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgesvd_Tile_Async( jobu, jobvt, A, S, T, U, LDU, VT, LDVT, sequence, &request );

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
 *  CHAMELEON_zgesvd_Tile_Async - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Non-blocking equivalent of CHAMELEON_zgesvd_Tile().
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
 * @sa CHAMELEON_zgesvd
 * @sa CHAMELEON_zgesvd_Tile
 * @sa CHAMELEON_cgesvd_Tile_Async
 * @sa CHAMELEON_dgesvd_Tile_Async
 * @sa CHAMELEON_sgesvd_Tile_Async
 *
 */
int CHAMELEON_zgesvd_Tile_Async( cham_job_t jobu, cham_job_t jobvt,
                                 CHAM_desc_t *A,
                                 double *S,
                                 CHAM_desc_t *T,
                                 CHAMELEON_Complex64_t *U, int LDU,
                                 CHAMELEON_Complex64_t *VT, int LDVT,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_desc_t descA;
    CHAM_desc_t descT;
    CHAM_desc_t D, *Dptr = NULL;
    double *E;
    int M, N, MINMN;

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();

    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgesvd_Tile_Async", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    assert( chamctxt->scheduler != RUNTIME_SCHED_PARSEC );
    if ( sequence == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgesvd_Tile_Async", "NULL sequence" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgesvd_Tile_Async", "NULL request" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if ( sequence->status == CHAMELEON_SUCCESS ) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED );
    }

    /* Check descriptors for correctness */
    if ( chameleon_desc_check(A) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "invalid first descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    } else {
        descA = *A;
    }
    if ( chameleon_desc_check(T) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "invalid fourth descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if ( (jobu != ChamNoVec) && (jobu != ChamAllVec) ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "illegal value of jobu" );
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if ( (jobvt != ChamNoVec) && (jobvt != ChamAllVec) ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "illegal value of jobvt" );
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    if ( descA.nb != descA.mb ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "only square tiles supported" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }

    M     = descA.m;
    N     = descA.n;
    MINMN = chameleon_min(M, N);
#if defined(CHAMELEON_COPY_DIAG)
    {
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &D, A->m, A->n );
        Dptr = &D;
    }
#endif

    E = malloc( sizeof(double) * MINMN );
    if ( E == NULL ) {
        chameleon_error( "CHAMELEON_zgesvd_Tile_Async", "malloc(E) failed" );
        free( E );
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }
    memset( E, 0, sizeof(double) * MINMN );


    /* Reduction to band + bidiagonal */
    chameleon_pzgebrd( 1, jobu, jobvt, &descA, &descT, Dptr,
                       U, LDU, VT, LDVT, E, S, sequence, request );

    /* Solve the bidiagonal SVD problem */
    /* On exit, U and VT are updated with bidiagonal matrix singular vectors */
#if !defined(CHAMELEON_SIMULATION)
    {
        cham_uplo_t uplo;
        int info, nru, ncvt;

        switch ( jobu ) {
            /* case ChamOVec : */
            case ChamAllVec :
                nru = M;
                break;
            /* case ChamSVec : */
            /*     nru = MINMN; */
            /*     break; */
            case ChamNoVec :
            default:
                nru = 0;
        }
        switch ( jobvt ) {
            /* case ChamOVec : */
            case ChamAllVec :
                ncvt = N;
                break;
            /* case ChamSVec : */
            /*     ncvt = MINMN; */
            /*     break; */
            case ChamNoVec :
            default:
                ncvt = 0;
        }
        uplo = M >= N ? ChamUpper : ChamLower;
        info = LAPACKE_zbdsqr( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), MINMN,
                               ncvt, nru, 0, S, E, VT, LDVT, U, LDU, NULL, 1 );
        if ( info != 0 ) {
            fprintf(stderr, "CHAMELEON_zgesvd_Tile_Async: LAPACKE_zbdsqr = %d\n", info );
        }
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    free( E );
    if ( Dptr != NULL ) {
        chameleon_desc_destroy( Dptr );
    }
    (void)D;
    return CHAMELEON_SUCCESS;
}
