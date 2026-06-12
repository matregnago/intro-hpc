/**
 *
 * @file testing_zcheck_polar_decomp.c
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
 * @date 2024-02-18
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
 * @brief Extends the check_zgeqrf and check_zortho to check the specific QR
 * factorization and Q generation used in Polar decompositions.
 *
 * [ A1 ] = [ Q1 ] [ AF1 ]
 * [ A2 ]   [ Q2 ]
 *
 *******************************************************************************
 *
 * @param[in] descA1
 *          The descriptor of the top of the matrix A.
 *
 * @param[in] descA2
 *          The descriptor of the bottom of the matrix A.
 *
 * @param[in] descQ1
 *          The descriptor of the top of the matrix Q generated from the
 *          QR factorization of the matrix A.
 *
 * @param[in] descQ2
 *          The descriptor of the bottom of the matrix Q generated from the
 *          QR factorization of the matrix A.
 *
 * @param[in] descAF1
 *          The descriptor of the top of the QR factorization of the matrix A that holds R.
 *
 * @retval 0  on success
 * @retval >0 on failure
 *
 *******************************************************************************
 */
int check_zgepdf_qr( run_arg_list_t *args, CHAM_desc_t *descA1, CHAM_desc_t *descA2, CHAM_desc_t *descQ1,
                     CHAM_desc_t *descQ2, CHAM_desc_t *descAF1 )
{
    int info_local, info_global;
    int M = (descQ1->m + descQ2->m);
    int N = descQ1->n;
    int K = descAF1->n;
    double result, Anorm, A1norm, A2norm, Rnorm;
    cham_fixdbl_t eps = testing_getaccuracy();
    CHAM_desc_t *descR, *subR, *subAF;

    /*
     * We exploit the fact that the lower matrix is supposed to be smaller than
     * the upper one, and the fact that K is always equal to N in this specific
     * problem.
     */
    assert( descAF1->m >= N );
    assert( K == N );

    /*
     * First check: || R - Q' A ||
     */
    descR = CHAMELEON_Desc_Copy( descAF1, CHAMELEON_MAT_ALLOC_TILE );

    /* Copy R(1:n,1:k) */
    subR  = chameleon_desc_submatrix( descR,   0, 0, N, K );
    subAF = chameleon_desc_submatrix( descAF1, 0, 0, N, K );

    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., subR );
    CHAMELEON_zlacpy_Tile( ChamUpper, subAF, subR );
    free( subAF );

    /*
     * Compute R(1:n,1:k) - Q(:,1:n)' * A(:,1:k)
     *       = R(1:n,1:k) - Q1(:,1:n)' * A1(:,1:k) - Q2(,1:n)' * A2(,1:k)
     */
    CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ1, descA1, 1., subR );
    CHAMELEON_zgemm_Tile( ChamConjTrans, ChamNoTrans, -1., descQ2, descA2, 1., subR );

    Rnorm  = CHAMELEON_zlange_Tile( ChamOneNorm, subR );
    A1norm = CHAMELEON_zlange_Tile( ChamOneNorm, descA1 );
    A2norm = CHAMELEON_zlange_Tile( ChamOneNorm, descA2 );
    Anorm  = A1norm + A2norm; /* This is an upper bound exact when A2 is diagonal due to OneNorm */
    result = Rnorm / ( (double)M * Anorm * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    /*
     * Second check: || I - Q' Q ||
     */
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., subR );

    /* Performs I - Q'Q = I - Q1'Q1 - Q2'Q2 */
    CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ1, 1., subR );
    CHAMELEON_zherk_Tile( ChamUpper, ChamConjTrans, -1., descQ2, 1., subR );

    /* Verifies the residual's norm */
    Rnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, ChamUpper, subR );
    result = Rnorm / ( (double)K * eps );

    run_arg_add_double( args, "||I-QQ'||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local++;
    }

    free( subR );
    CHAMELEON_Desc_Destroy( &descR );

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
 * @brief Checks if a Chameleon Polar Decomposition is correct.
 *
 * @Warning Check only the general case for now.
 *
 *******************************************************************************
 *
 * @param[in,out] descA
 *          The descriptor of the matrix A, on exit the matrix is modified.
 *
 * @param[in] descU
 *          The descriptor of the orthogonal polar factor of the decomposition.
 *
 * @param[in] descH
 *          The descriptor of the symmetric/hermitian polar factor of the decomposition.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxpd( run_arg_list_t *args,
                 CHAM_desc_t *descA, CHAM_desc_t *descU, CHAM_desc_t *descH )
{
    int info_local, info_global;
    double Anorm, Rnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();

    /* Compute ||A|| */
    Anorm = CHAMELEON_zlange_Tile( ChamFrobeniusNorm, descA );

    /* R = A - U * H */
    CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, 1., descU, descH, -1., descA );

    /* Compute ||R|| */
    Rnorm = CHAMELEON_zlange_Tile( ChamFrobeniusNorm, descA );

    result = Rnorm / (Anorm * eps);
    run_arg_add_double( args, "||A||",         Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else{
        info_local = 0;
    }

    /* Broadcasts the result from the main processus */
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
 * @brief Checks if a Polar Decomposition is correct.
 *
 * @Warning Check only the general case for now.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrices A and U.
 *
 * @param[in] N
 *          The number of columns of the matrices A and U and the order of the matrix H.
 *
 * @param[in,out] A
 *          The matrix A, on exit the matrix is modified.
 *
 * @param[in] U
 *          The matrix with the orthogonal polar factor of the decomposition.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and U.
 *
 * @param[in] H
 *          The matrix of the symmetric/hermitian polar factor of the decomposition.
 *
 * @param[in] LDH
 *          The leading dimension of the matrix H.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxpd_std( run_arg_list_t *args, int M, int N, CHAMELEON_Complex64_t *A,
                     CHAMELEON_Complex64_t *U, int LDA, CHAMELEON_Complex64_t *H, int LDH )
{
    int          info;
    CHAM_desc_t *descA, *descU, *descH;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descU, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descH, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDH, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, U, LDA, descU );
    CHAMELEON_zLap2Desc( ChamUpperLower, H, LDH, descH );

    info = check_zxxpd( args, descA, descU, descH );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descU );
    CHAMELEON_Desc_Destroy( &descH );

    return info;
}

#endif /* defined(CHAMELEON_SIMULATION) */
