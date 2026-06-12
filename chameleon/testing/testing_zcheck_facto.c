/**
 *
 * @file testing_zcheck_facto.c
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

#ifndef max
#define max( _a_, _b_ ) ( (_a_) > (_b_) ? (_a_) : (_b_) )
#endif

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed product U*U' or L'*L result with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether the upper or lower triangle of A is stored.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descAAt
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zlauum( run_arg_list_t *args, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descAAt )
{
    int info_local, info_global;
    int N       = descA->n;
    cham_fixdbl_t eps = testing_getaccuracy();
    double result, Anorm, AAtnorm, Rnorm;
    CHAM_desc_t *descAt;

    Anorm   = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, ChamNonUnit, descA );
    AAtnorm = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, ChamNonUnit, descAAt );

    if ( uplo == ChamUpper ) {
        descAt = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zlaset_Tile( ChamLower, 0., 0., descAt );
        CHAMELEON_zlacpy_Tile( ChamUpper, descA, descAt );

        /* Computes U * U' */
        CHAMELEON_ztrmm_Tile( ChamRight, ChamUpper, ChamConjTrans, ChamNonUnit, 1., descA, descAt );
    }
    else {
        descAt = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zlaset_Tile( ChamUpper, 0., 0., descAt );
        CHAMELEON_zlacpy_Tile( ChamLower, descA, descAt );

        /* Computes L' * L */
        CHAMELEON_ztrmm_Tile( ChamLeft, ChamLower, ChamConjTrans, ChamNonUnit, 1., descA, descAt );
    }

    /* Computes AAt - A * A' */
    CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, -1., descAAt, 1., descAt );

    Rnorm = CHAMELEON_zlantr_Tile( ChamMaxNorm, uplo, ChamNonUnit, descAt );

    CHAMELEON_Desc_Destroy( &descAt );

    /* Compares the residual's norm */
    result = Rnorm / ( Anorm * Anorm * N * eps );
    if (  isnan(AAtnorm) || isinf(AAtnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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
 * @brief Compares two product U*U' or L'*L result.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Whether the upper or lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrices A and AAt.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] AAt
 *          The matrix of the computed result.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and AAt.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zlauum_std( run_arg_list_t *args, cham_uplo_t uplo, int N, CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *AAt, int LDA )
{
    int          info;
    CHAM_desc_t *descA;
    CHAM_desc_t *descAAt;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,   CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAAt, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( uplo, A,   LDA, descA   );
    CHAMELEON_zLap2Desc( uplo, AAt, LDA, descAAt );

    info = check_zlauum( args, uplo, descA, descAAt );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descAAt );

    (void)args;
    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if a Chameleon computed factorization is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descLU
 *          The descriptor of the matrix of the Chameleon factorisation of the matrix A.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxtrf( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, CHAM_desc_t *descA, CHAM_desc_t *descLU )
{
    int info_local, info_global;
    int M      = descA->m;
    int N      = descA->n;
    double Anorm, Rnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();

    CHAM_desc_t *descL, *descU;
    cham_trans_t transL = ChamNoTrans;
    cham_trans_t transU = ChamNoTrans;

    descL = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
    descU = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descL );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 0., descU );

    switch ( uplo ) {
    case ChamUpper:
#if defined(PRECISION_z) || defined(PRECISION_c)
        transL = ( matrix_type == ChamHermitian ) ? ChamConjTrans : ChamTrans;
#else
        transL = ChamTrans;
#endif
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descL );
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descU );
        break;
    case ChamLower:
#if defined(PRECISION_z) || defined(PRECISION_c)
        transU = ( matrix_type == ChamHermitian ) ? ChamConjTrans : ChamTrans;
#else
        transU = ChamTrans;
#endif
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descL );
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descU );
        break;
    case ChamUpperLower:
    default:
        CHAMELEON_zlacpy_Tile( ChamLower, descLU, descL );
        CHAMELEON_zlaset_Tile( ChamUpper, 0., 1., descL );
        CHAMELEON_zlacpy_Tile( ChamUpper, descLU, descU );
    }

    switch ( matrix_type ) {
    case ChamGeneral: {
        CHAM_desc_t *subL, *subU;
        subL = chameleon_desc_submatrix( descL, 0, 0, M, chameleon_min(M, N) );
        subU = chameleon_desc_submatrix( descU, 0, 0, chameleon_min(M, N), N );

        Anorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., subL, subU, 1., descA );
        Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descA );

        free( subL );
        free( subU );
    }
        break;

        /* WARNING: A must be fully initialized and not just the correct triangular part */
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        Anorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., descL, descU, 1., descA );
        Rnorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA );
        break;
#endif

    case ChamSymmetric:
        Anorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA );
        CHAMELEON_zgemm_Tile( transL, transU, -1., descL, descU, 1., descA );
        Rnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA );
        break;

    default:
        fprintf(stderr, "check_zxxtrf: matrix_type(%d) unsupported\n", matrix_type );
        return 1;
    }

    result = Rnorm / ( Anorm * N * eps );
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||A-fact(A)||", Rnorm );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else{
        info_local = 0;
    }

    CHAMELEON_Desc_Destroy( &descL );
    CHAMELEON_Desc_Destroy( &descU );

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
 * @brief Checks if a core function computed factorization is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] M
 *          The number of rows of the matrices A and LU.
 *
 * @param[in] N
 *          The number of columns of the matrices A and LU.
 *
 * @param[in] A
 *          The symmetric matrix A.
 *
 * @param[in] LU
 *          The matrix with the computed factorisation of the matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and LU.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zxxtrf_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, int M, int N,
                      CHAMELEON_Complex64_t *A, CHAMELEON_Complex64_t *LU, int LDA )
{
    int          info;
    CHAM_desc_t *descA;
    CHAM_desc_t *descLU;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA,  CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descLU, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, 1, 1 );

    CHAMELEON_zLap2Desc( ChamUpperLower, A,  LDA, descA  );
    CHAMELEON_zLap2Desc( ChamUpperLower, LU, LDA, descLU );

    info = check_zxxtrf( args, matrix_type, uplo, descA, descLU );

    CHAMELEON_Desc_Destroy( &descA  );
    CHAMELEON_Desc_Destroy( &descLU );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the  linear solution of op(A) * x = b is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] trans
 *          Whether the A matrix is non transposed, tranposed or conjugate transposed.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descX
 *          The descriptor of the matrix X.
 *
 * @param[inout] descB
 *          The descriptor of the matrix B = A*X. On exit, it contains the remainder from A*x-B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsolve( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_trans_t trans, cham_uplo_t uplo,
                  CHAM_desc_t *descA, CHAM_desc_t *descX, CHAM_desc_t *descB )
{
    int info_local, info_global;
    int M                = descA->m;
    int N                = descA->n;
    double Anorm, Bnorm, Xnorm, Rnorm, result = 0;
    cham_fixdbl_t eps = testing_getaccuracy();
    cham_normtype_t norm = (trans == ChamNoTrans) ? ChamOneNorm : ChamInfNorm;

    /* Computes the norms */
    Bnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    Xnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descX );

    switch ( matrix_type ) {
    case ChamGeneral:
        Anorm = CHAMELEON_zlange_Tile( norm, descA );
        CHAMELEON_zgemm_Tile( trans, ChamNoTrans, -1., descA, descX, 1., descB );
        break;

#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        Anorm = CHAMELEON_zlanhe_Tile( norm, uplo, descA );
        CHAMELEON_zhemm_Tile( ChamLeft, uplo, -1., descA, descX, 1., descB );
        break;
#endif

    case ChamSymmetric:
        Anorm = CHAMELEON_zlansy_Tile( norm, uplo, descA );
        CHAMELEON_zsymm_Tile( ChamLeft, uplo, -1., descA, descX, 1., descB );
        break;

    default:
        fprintf(stderr, "check_zsolve: matrix_type(%d) unsupported\n", matrix_type );
        return 1;
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descB );
    result = Rnorm / ( Anorm * Xnorm * max( M, N ) * eps );

    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||X||", Xnorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||Ax-b||", Rnorm );
    run_arg_add_double( args, "||Ax-b||/N/eps/(||A||||x||+||b||", result );

    if ( isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
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

    return info_global;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the  linear solution of op(A) * x = b is correct.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] trans
 *          Whether the A matrix is non transposed, tranposed or conjugate transposed.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] N
 *          The order of the matrix A and the number of rows of the matrices X and B.
 *
 * @param[in] NRHS
 *          The number of columns of the matrices X and B.
 *
 * @param[in] A
 *          The symmetric matrix A.
 *
 * @param[in] LDA
 *          The leading dimenson of the matrix A.
 *
 * @param[in] X
 *          The matrix X.
 *
 * @param[inout] B
 *          The matrix B = A*X. On exit, it contains the remainder from A*x-B.
 *
 * @param[in] LDB
 *          The leading dimension of the matrices X and B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsolve_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_trans_t trans, cham_uplo_t uplo, int N, int NRHS,
                      CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *X, CHAMELEON_Complex64_t *B, int LDB )
{
    int          info;
    CHAM_desc_t *descA, *descX, *descB;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N,    0, 0, N, N,    1, 1 );
    CHAMELEON_Desc_Create(
        &descB, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, 1, 1 );
    CHAMELEON_Desc_Create(
        &descX, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, N, NRHS, 1, 1 );

    CHAMELEON_zLap2Desc( uplo,           A, LDA, descA );
    CHAMELEON_zLap2Desc( ChamUpperLower, B, LDB, descB );
    CHAMELEON_zLap2Desc( ChamUpperLower, X, LDB, descX );

    info = check_zsolve( args, matrix_type, trans, uplo, descA, descX, descB );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descX );

    return info;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the matrix A0 is the inverse of Ai.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether they are upper triangular matrices or lower triangular matrices.
 *
 * @param[in] diag
 *          Whether they are unitary diagonal matrices or not.
 *
 * @param[in] descA0
 *          The descriptor of the matrix A0.
 *
 * @param[in] descAi
 *          The descriptor of the matrix Ai.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrtri( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                  CHAM_desc_t *descA0, CHAM_desc_t *descAi )
{
    int          info_local, info_global;
    cham_uplo_t  uplo_inv;
    CHAM_desc_t *descI, *descB = NULL;
    double       Rnorm, Anorm, Ainvnorm, result;
    cham_fixdbl_t eps = testing_getaccuracy();
    int          N   = descA0->m;

    /* Creates an identity matrix */
    descI = CHAMELEON_Desc_Copy( descA0, CHAMELEON_MAT_ALLOC_TILE );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., descI );

    /* Calculates the residual I - A*(A**-1) */
    switch ( matrix_type ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        /* Ainv comes from potri and is hermitian */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descA0 );
        Ainvnorm = CHAMELEON_zlanhe_Tile( ChamOneNorm, uplo, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_ztradd_Tile( uplo_inv, ChamConjTrans, 1., descAi, 0., descB );
        CHAMELEON_zlacpy_Tile( uplo, descAi, descB );

        CHAMELEON_zhemm_Tile( ChamLeft, uplo, -1., descA0, descB, 1., descI );
        break;
#endif

    case ChamSymmetric:
        /* Ainv comes from potri and is symmetric */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descA0 );
        Ainvnorm = CHAMELEON_zlansy_Tile( ChamOneNorm, uplo, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_ztradd_Tile( uplo_inv, ChamTrans, 1., descAi, 0., descB );
        CHAMELEON_zlacpy_Tile( uplo, descAi, descB );

        CHAMELEON_zsymm_Tile( ChamLeft, uplo, -1., descA0, descB, 1., descI );
        break;

    case ChamTriangular:
        /* Ainv comes from trtri */
        assert( uplo != ChamUpperLower );

        Anorm    = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, diag, descA0 );
        Ainvnorm = CHAMELEON_zlantr_Tile( ChamOneNorm, uplo, diag, descAi );

        /*
         * Expand Ainv into a full matrix and call ZHEMM to multiply
         * Ainv on the left by A.
         */
        uplo_inv = ( uplo == ChamUpper ) ? ChamLower : ChamUpper;
        descB = CHAMELEON_Desc_Copy( descAi, CHAMELEON_MAT_ALLOC_TILE );

        if ( diag == ChamUnit ) {
            //CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, 1., descAi, 0., descB );
            CHAMELEON_zlacpy_Tile( uplo, descAi, descB );
            CHAMELEON_zlaset_Tile( uplo_inv, 0., 1., descB );
        }
        else {
            CHAMELEON_zlaset_Tile( uplo_inv, 0., 1., descB );
            CHAMELEON_zlacpy_Tile( uplo, descAi, descB );
            //CHAMELEON_ztradd_Tile( uplo, ChamNoTrans, 1., descAi, 0., descB );
        }

        /* Computes - A * A^-1 */
        CHAMELEON_ztrmm_Tile( ChamLeft, uplo, ChamNoTrans, diag, -1., descA0, descB );
        /* Computes I - A * A^-1 */
        CHAMELEON_zgeadd_Tile( ChamNoTrans, 1., descB, 1., descI );
        break;

    case ChamGeneral:
    default:
        /* Ainv comes from getri */
        assert( uplo == ChamUpperLower );

        Anorm    = CHAMELEON_zlange_Tile( ChamOneNorm, descA0 );
        Ainvnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descAi );

        CHAMELEON_zgemm_Tile( ChamNoTrans, ChamNoTrans, -1., descA0, descAi, 1., descI );
        break;
    }

    Rnorm = CHAMELEON_zlange_Tile( ChamOneNorm, descI );

    /* Compares the residual's norm */
    result = Rnorm / ( Anorm * Ainvnorm * N * eps );
    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_local = 1;
    }
    else {
        info_local = 0;
    }

    CHAMELEON_Desc_Destroy( &descI );
    if ( descB != NULL ) {
        CHAMELEON_Desc_Destroy( &descB );
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
 * @brief Checks if the matrix A0 is the inverse of Ai.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether they are upper triangular matrices or lower triangular matrices.
 *
 * @param[in] diag
 *          Whether they are unitary diagonal matrices or not.
 *
 * @param[in] N
 *          the order of the matrices A0 and Ai.
 *
 * @param[in] A0
 *          The matrix A0.
 *
 * @param[in] Ai
 *          The matrix Ai.
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A0 and Ai.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrtri_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_diag_t diag,
                      int N, CHAMELEON_Complex64_t *A0, CHAMELEON_Complex64_t *Ai, int LDA )
{
    int          info;
    CHAM_desc_t *descA0, *descAi;

    int nb;
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );

    CHAMELEON_Desc_Create(
        &descA0, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );
    CHAMELEON_Desc_Create(
        &descAi, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, 1, 1 );

    CHAMELEON_zLap2Desc( uplo, A0, LDA, descA0 );
    CHAMELEON_zLap2Desc( uplo, Ai, LDA, descAi );

    info = check_ztrtri( args, matrix_type, uplo, diag, descA0, descAi );

    CHAMELEON_Desc_Destroy( &descA0 );
    CHAMELEON_Desc_Destroy( &descAi );

    return info;
}

#endif /* defined(CHAMELEON_SIMULATION) */
