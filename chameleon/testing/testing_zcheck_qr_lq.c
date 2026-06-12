/**
 *
 * @file testing_zcheck_qr_lq.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @author Alycia Lisito
 * @author Pierre Esterie
 * @date 2024-11-13
 * @precisions normal z -> c d s
 *
 */
#include "../control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <chameleon.h>

#if !defined(CHAMELEON_SIMULATION)

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a matrix is orthogonal.
 *
 *******************************************************************************
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zortho( run_arg_list_t *args, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M      = descQ->m;
    int N      = descQ->n;
    int minMN  = chameleon_min(M, N);
    double result, normR;
    cham_fixdbl_t eps = testing_getaccuracy();
    CHAM_desc_t *descI, *subI;

    /* Builds the identity */
    descI = CHAMELEON_Desc_Copy( descQ, CHAMELEON_MAT_ALLOC_TILE );
    subI = chameleon_desc_submatrix( descI, 0, 0, minMN, minMN );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., subI );

    /* Performs Id - Q'Q */
    if ( M >= N ) {
        CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ, 1., subI );
    }
    else {
        CHAMELEON_zherk_Tile( ChamUpper, ChamNoTrans, -1., descQ, 1., subI );
    }

    /* Verifies the residual's norm */
    normR = CHAMELEON_zlansy_Tile( ChamOneNorm, ChamUpper, subI );
    result = normR / ( (double)minMN * eps );

    run_arg_add_double( args, "||I-QQ'||", normR );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    free( subI );
    CHAMELEON_Desc_Destroy( &descI );

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a matrix is orthogonal.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix Q.
 *
 * @param[in] N
 *          The number of columns of the matrix Q.
 *
 * @param[in] descQ
 *          The matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zortho_std( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descQ, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, Q, LDQ, descQ );

    info = check_zortho( args, descQ );

    CHAMELEON_Desc_Destroy( &descQ );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqf( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M = descQ->m;
    int N = descQ->n;
    int K = chameleon_min( descA->m, descA->n );
    double result, Anorm, Rnorm;
    cham_fixdbl_t eps = testing_getaccuracy();
    CHAM_desc_t *descL;
    int full_lq = ( M == N ) ? 1 : 0;

    assert( descA->n == N );
    assert( descA->m == descAF->m );

    descL = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

    if ( full_lq ) {
        /*
         * Cas lapack zlqt01.f
         * A full N-by-N Q has been generated
         */
        assert( descAF->n == N );

        /* Copy L */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descL );
        CHAMELEON_zlacpy_Tile( ChamLower, descAF, descL );

        /* Compute L - A * Q' */
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamConjTrans, -1., descA, descQ, 1., descL );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descL );
    }
    else {
        /*
         * Cas lapack zlqt02.f
         * A partial Q has been generated (K < min(M, N))
         */
        CHAM_desc_t *subL, *subAF;

        assert( descAF->n >= M );
        assert( N >= M );

        /* Copy L(1:k,1:m) */
        subL  = chameleon_desc_submatrix( descL,  0, 0, K, M );
        subAF = chameleon_desc_submatrix( descAF, 0, 0, K, M );

        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subL );
        CHAMELEON_zlacpy_Tile( ChamLower, subAF, subL );

        /* Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)' */
        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamConjTrans, -1., descA, descQ, 1., subL );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, subL );

        free( subL );
        free( subAF );
    }

    CHAMELEON_Desc_Destroy(&descL);

    Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    result = Rnorm / ( (double)N * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] K
 *          The number of elementary reflectors whose product defines the matrix Q.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] Q
 *          The matrix Q generated with a call to ungqr and the computed factorisation of the matrix A (descAF).
 *
 * @param[in] LDQ
 *          The leading dimension fo the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqf_std( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descA, *descAF, *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAF, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, AF, LDA, descAF );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ  );

    info = check_zgelqf( args, descA, descAF, descQ );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAF );
    CHAMELEON_Desc_Destroy( &descQ  );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrf( run_arg_list_t *args, CHAM_desc_t *descA, CHAM_desc_t *descAF, CHAM_desc_t *descQ )
{
    int info_local, info_global;
    int M = descQ->m;
    int N = descQ->n;
    int K = chameleon_min( descA->m, descA->n );
    double result, Anorm, Rnorm;
    cham_fixdbl_t eps = testing_getaccuracy();
    CHAM_desc_t *descR;
    int full_qr = ( M == N ) ? 1 : 0;

    assert( descA->m == M );
    assert( descA->n == descAF->n );

    descR = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

    if ( full_qr ) {
        /*
         * Cas lapack zqrt01.f
         * A full M-by-M Q has been generated
         */
        assert( descAF->m == M );

        /* Copy R */
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descR );
        CHAMELEON_zlacpy_Tile( ChamUpper, descAF, descR );

        /* Compute R - Q'*A */
        CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ, descA, 1., descR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descR );
    }
    else {
        /*
         * Cas lapack zqrt02.f
         * A partial Q has been generated (K < min(M, N))
         */
        CHAM_desc_t *subR, *subAF;

        assert( descAF->m >= N );
        assert( N <= M );

        /* Copy R(1:n,1:k) */
        subR  = chameleon_desc_submatrix( descR,  0, 0, N, K );
        subAF = chameleon_desc_submatrix( descAF, 0, 0, N, K );

        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subR );
        CHAMELEON_zlacpy_Tile( ChamUpper, subAF, subR );

        /* Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k) */
        CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ, descA, 1., subR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, subR );

        free( subR );
        free( subAF );
    }

    CHAMELEON_Desc_Destroy(&descR);

    Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    result = Rnorm / ( (double)M * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Reduces the result on all processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a linear solution is correct.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] K
 *          The number of elementary reflectors whose product defines the matrix Q.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] descAF
 *          The descriptor of the matrix of the Chameleon computed factorisation of the matrix A.
 *
 * @param[in] Q
 *          The matrix Q generated with a call to ungqr and
 *          the computed factorisation of the matrix A (descAF).
 *
 * @param[in] LDQ
 *          The leading dimension fo the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrf_std( run_arg_list_t *args, int M, int N, int K, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AF, int LDA, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int          info;
    CHAM_desc_t *descA, *descAF, *descQ;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAF, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, 1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, M, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, AF, LDA, descAF );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ  );

    info = check_zgeqrf( args, descA, descAF, descQ );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAF );
    CHAMELEON_Desc_Destroy( &descQ  );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the decomposition is correct.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Whether the matrix Q appears on the left or on the right of the product.
 *
 * @param[in] trans
 *          Whether the matrix Q is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descC
 *          The descriptor of the matrix C.
 *
 * @param[in] descQ
 *          The descriptor of the matrix Q.
 *
 * @param[in] descCC
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zqc( run_arg_list_t *args, cham_side_t side, cham_trans_t trans,
               CHAM_desc_t *descC, CHAM_desc_t *descQ, CHAM_desc_t *descCC )
{
    int info_local, info_global;
    int M = descQ->m;
    double Cnorm, Qnorm, CCnorm, Rnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();

    Cnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, descC );
    Qnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, descQ );
    CCnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descCC );

    if ( side == ChamLeft ) {
        CHAMELEON_zgemm_Tile( trans, ChamNoTrans, -1., descQ, descC, 1., descCC );
    }
    else {
        CHAMELEON_zgemm_Tile( ChamNoTrans, trans, -1., descC, descQ, 1., descCC );
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descCC );
    result = Rnorm / ( M * Cnorm * eps );

    run_arg_add_double( args, "||A||", Qnorm );
    run_arg_add_double( args, "||B||", CCnorm );
    run_arg_add_double( args, "||C||", Cnorm );
    run_arg_add_double( args, "||R||", Rnorm );

    if ( isnan(CCnorm) || isinf(CCnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)Qnorm;
    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the decomposition is correct.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Whether the matrix Q appears on the left or on the right of the product.
 *
 * @param[in] trans
 *          Whether the matrix Q is transposed, conjugate transposed or not transposed.
 *
 * @param[in] M
 *          The number of rows of the matrices C and CC and the order of the matrix Q if side = ChamLeft.
 *
 * @param[in] N
 *          The number of columns of the matrices C and CC and the order of the matrix Q if side = ChamRight.
 *
 * @param[in] C
 *          The matrix C.
 *
 * @param[in] CC
 *          The matrix of the computed result.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and CC.
 *
 * @param[in] Q
 *          The matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the matrix Q.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zqc_std( run_arg_list_t *args, cham_side_t side, cham_trans_t trans, int M, int N,
                   CHAMELEON_Complex64_t *C, CHAMELEON_Complex64_t *CC, int LDC, CHAMELEON_Complex64_t *Q, int LDQ )
{
    int info;
    int Qm = ( side == ChamLeft ) ? M : N;
    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAM_desc_t *descC, *descCC, *descQ;

    CHAMELEON_Desc_Create(
        &descC,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N,  0, 0, M,  N,  1, 1 );
    CHAMELEON_Desc_Create(
        &descCC, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N,  0, 0, M,  N,  1, 1 );
    CHAMELEON_Desc_Create(
        &descQ,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDQ, Qm, 0, 0, Qm, Qm, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, C,  LDC, descC );
    CHAMELEON_zLap2Desc( ChamUpperLower, CC, LDC, descCC );
    CHAMELEON_zLap2Desc( ChamUpperLower, Q,  LDQ, descQ );

    info = check_zqc( args, side, trans, descC, descQ, descCC );

    CHAMELEON_Desc_Destroy( &descC );
    CHAMELEON_Desc_Destroy( &descCC );
    CHAMELEON_Desc_Destroy( &descQ );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] Bnorm
 *          TODO
 *
 * @param[in] descR
 *          The descriptor of the matrix R.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgeqrs( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR )
{
    int info_local, info_global, nb;
    int M = descA->m;
    int N = descA->n;
    int NRHS = descR->n;
    int maxMNK = chameleon_max( M, chameleon_max( N, NRHS ) );
    double Rnorm, result;
    double Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    cham_fixdbl_t eps = testing_getaccuracy();

    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    if ( trans == ChamNoTrans ) {
        CHAM_desc_t *descRR;
        /*
         * Corresponds to lapack/testings/lin/[sdcz]qrt17.f
         *
         * ZQRT17 computes the ratio
         *
         *    || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
         *
         * where R = op(A)*X - B, op(A) is A or A', and alpha = ||B||
         *
         */
        CHAMELEON_Desc_Create( &descRR, NULL, ChamComplexDouble, nb, nb, nb*nb,
                               NRHS, N, 0, 0, NRHS, N, chameleon_desc_datadist_get_iparam(descA, 0), chameleon_desc_datadist_get_iparam(descA, 1) );

        CHAMELEON_zgemm_Tile( ChamConjTrans, trans, 1., descR, descA, 0., descRR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descRR );
        result = Rnorm / (Anorm * Bnorm * eps * maxMNK);
        CHAMELEON_Desc_Destroy( &descRR );
    }
    else {
        /*
         * To implement this test, we need to look at LAPACK working note 41, page 29
         * and more especially to lapack/testings/lin/[sdcz]qrt14.f
         */
        fprintf(stderr, "GEQRS testing not implemented with M >= N when transA = ChamConjTrans\n");
        return 0;
    }

    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] Bnorm
 *          TODO
 *
 * @param[in] descR
 *          The descriptor of the matrix R.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgelqs( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, double Bnorm, CHAM_desc_t *descR )
{
    int info_local, info_global, nb;
    int M = descA->m;
    int N = descA->n;
    int NRHS = descR->n;
    int maxMNK = chameleon_max( M, chameleon_max( N, NRHS ) );
    double Rnorm, result;
    double Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
    cham_fixdbl_t eps = testing_getaccuracy();

    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    if ( trans == ChamNoTrans ) {
        /*
         * To implement this test, we need to look at LAPACK working note 41, page 29
         * and more especially to lapack/testings/lin/[sdcz]lqt14.f
         */
        fprintf(stderr, "GELQS testing not implemented with N > M when transA = ChamNoTrans\n");
        return 0;
    }
    else {
        CHAM_desc_t *descRR;
        /*
         * Corresponds to lapack/testings/lin/[sdcz]qrt17.f
         *
         * ZQRT17 computes the ratio
         *
         *    || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
         *
         * where R = op(A)*X - B, op(A) is A or A', and alpha = ||B||
         *
         */
        CHAMELEON_Desc_Create( &descRR, NULL, ChamComplexDouble, nb, nb, nb*nb, NRHS, M, 0, 0, NRHS, M,
                               chameleon_desc_datadist_get_iparam(descA, 0), chameleon_desc_datadist_get_iparam(descA, 1) );

        CHAMELEON_zgemm_Tile( ChamConjTrans, trans, 1., descR, descA, 0., descRR );

        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descRR );
        result = Rnorm / (Anorm * Bnorm * eps * maxMNK);
        CHAMELEON_Desc_Destroy( &descRR );
    }

    if ( isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Allreduce( &info_local, &info_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#else
    info_global = info_local;
#endif

    (void)args;
    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descX
 *          The descriptor of the matrix X.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgels( run_arg_list_t *args, cham_trans_t trans, CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB )
{
    int info_solution;
    int M = descA->m;
    int N = descA->n;
    double Bnorm = CHAMELEON_zlange_Tile( ChamInfNorm, descB );

    info_solution = check_zsolve( args, ChamGeneral, trans, ChamUpperLower,
                                  descA, descX, descB );

    if ( M >= N ) {
        info_solution = check_zgeqrs( args, trans, descA, Bnorm, descB );
    }
    else {
        info_solution = check_zgelqs( args, trans, descA, Bnorm, descB );
    }

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Whether the matrix A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] X
 *          The matrix X.
 *
 * @param[in] LDX
 *          The leading dimension of the matrix X.
 *
 * @param[in] B
 *          The matrix B.
 *
 * @param[in] LDB
 *          The leading dimension fo the matrix B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgels_std( run_arg_list_t *args, cham_trans_t trans, int M, int N, int NRHS, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, int LDX, CHAMELEON_Complex64_t *B, int LDB )
{
    int info_solution;
    CHAM_desc_t *descA, *descX, *descB;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );
    int Xm = (trans == ChamNoTrans) ? N : M;
    int Bm = (trans == ChamNoTrans) ? M : N;

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N,    0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descX, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDX, NRHS, 0, 0, Xm, NRHS, 1, 1 );
    CHAMELEON_Desc_Create(
        &descB, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, Bm, NRHS, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, X, LDX, descX );
    CHAMELEON_zLap2Desc( ChamUpperLower, B, LDB, descB );

    info_solution = check_zgels( args, trans, descA, descX, descB );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descB );

    return info_solution;
}

#endif /* defined(CHAMELEON_SIMULATION) */
