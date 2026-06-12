/**
 *
 * @file testing_zcheck_blas.c
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
 * @date 2025-01-29
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

#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas.h>
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
 * @brief Compares two core function computed products.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] transB
 *          Whether the second product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] M
 *          The number of rows of the matrices A and C.
 *
 * @param[in] N
 *          The number of columns of the matrices B and C.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the  umber of rows of the matrxi B.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] B
 *          The matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The matrix Cref.
 *
 * @param[in] C
 *          The matrix of the computed result alpha*A*B+beta*C.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and Cref.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgemm_std( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, int M, int N, int K, CHAMELEON_Complex64_t *A, int LDA,
                     CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC )
{
    int info_solution           = 0;
    double Anorm, Bnorm, Crefnorm, Rnorm, result;
    CHAMELEON_Complex64_t mzone = -1.0;

    /* Calculates the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', M, K, A, LDA );
    } else {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, M, A, LDA );
    }
    if ( transB == ChamNoTrans ) {
        Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', K, N, B, LDB );
    } else {
        Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', N, K, B, LDB );
    }

    /* Computes the norms for comparing */
    Crefnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    cham_fixdbl_t eps = testing_getaccuracy();

    /* Makes the multiplication with the core function */
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                 CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
    cblas_zaxpy( LDC * N, CBLAS_SADDR(mzone), C, 1, Cref, 1 );

    /* Calculates the norm with the core function's result */
    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    if ( ( alpha != 0. ) || (beta != 0. ) ) {
        result = Rnorm / ( ( cabs(alpha) * max(Anorm, Bnorm) + cabs(beta) * Crefnorm ) * K * eps );
    }
    else {
        result = Rnorm;
    }
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||C||", Crefnorm );
    run_arg_add_double( args, "||R||", Rnorm );

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] transB
 *          Whether the second product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix Cref.
 *
 * @param[in] descC
 *          The descriptor of the matrix of the Chameleon computed result alpha*A*B+beta*C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgemm( run_arg_list_t *args, cham_trans_t transA, cham_trans_t transB, CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA,
                 CHAM_desc_t *descB, CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution           = 0;
    int M                       = descC->m;
    int N                       = descC->n;
    int K                       = (transA != ChamNoTrans)? descA->m : descA->n;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *C    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    int An                      = ( transA == ChamNoTrans ) ? K : M;
    int LDA                     = ( transA == ChamNoTrans ) ? M : K;
    int Bn                      = ( transB == ChamNoTrans ) ? N : K;
    int LDB                     = ( transB == ChamNoTrans ) ? K : N;
    int LDC                     = M;

    /* Creates the LAPACK version of the matrices */
    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDA * An );
        B    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDB * Bn );
        Cref = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDC * N  );
        C    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDC * N  );
    }

    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {

        info_solution = check_zgemm_std( args, transA, transB, alpha, M, N, K, A, LDA, B, LDB, beta, Cref, C, LDC );

        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed symmetric products.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] side
 *          Whether the symmetric matrix A appears on the left or right in the operation.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] M
 *          The order of the matrix A and number of rows of the matrices B and C.
 *
 * @param[in] N
 *          The number of columns of the matrices B and C.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The symmetric matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] B
 *          The matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The matrix Cref.
 *
 * @param[in] C
 *          The matrix of the computed result alpha*A*B+beta*C or alpha*B*A+beta*C.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices C and Cref.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsymm_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo, int M, int N,
                     CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB,
                     CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref, CHAMELEON_Complex64_t *C, int LDC )
{
    int info_solution = 0;
    int An;
    double Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    CHAMELEON_Complex64_t mzone = -1.0;
    double *work;
    char normTypeA, normTypeB;

    if ( side == ChamLeft ) {
        normTypeA = 'I';
        normTypeB = 'O';
        An = M;
    }
    else {
        normTypeA = 'O';
        normTypeB = 'I';
        An = N;
    }

    work = malloc( sizeof(double) * An );
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Anorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, normTypeA, chameleon_lapack_const(uplo), An, A, LDA, work );
    }
    else
#endif
    {
        Anorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, normTypeA, chameleon_lapack_const(uplo), An, A, LDA, work );
    }
    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, normTypeB, M, N, B, LDB );
    free( work );

    /* Computes the norms for comparing */
    Crefnorm  = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );
    Cchamnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, C,    LDC, NULL );

    cham_fixdbl_t eps = testing_getaccuracy();

    /* Makes the multiplication with the core function */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        cblas_zhemm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                        M, N, CBLAS_SADDR(alpha),
                        A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
    }
    else
#endif
    {
        cblas_zsymm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB,
                     CBLAS_SADDR(beta), Cref, LDC );
    }
    cblas_zaxpy(LDC * N, CBLAS_SADDR(mzone), C, 1, Cref, 1);

    /* Computes the norm with the core function's result */
    Clapacknorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );
    Rnorm       = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, Cref, LDC, NULL );

    if ( ( alpha != 0. ) || (beta != 0. ) ) {
        result = Rnorm / ((cabs(alpha) * max(Anorm, Bnorm) + cabs(beta) * Crefnorm) * An * eps);
    }
    else {
        result = Rnorm;
    }

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    (void)Clapacknorm;
    (void)Cchamnorm;
    (void)matrix_type;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed symmetric product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] side
 *          Whether the symmetric matrix A appears on the left or right in the operation.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the symmetric matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the matrix Cref.
 *
 * @param[in] descC
 *          The descriptor of the matrix of the Chameleon computed result alpha*A*B+beta*C or alpha*B*A+beta*C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsymm( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_side_t side, cham_uplo_t uplo,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution = 0;
    int M             = descC->m;
    int N             = descC->n;
    int rank          = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;
    int An                      = ( side == ChamLeft )? M : N;
    int LDA                     = An;
    int LDB                     = M;
    int LDC                     = M;

    if ( rank == 0 ) {
        A    = ( CHAMELEON_Complex64_t * ) malloc( sizeof( CHAMELEON_Complex64_t ) * LDA * An );
        B    = ( CHAMELEON_Complex64_t * ) malloc( sizeof( CHAMELEON_Complex64_t ) * LDB * N  );
        Cref = ( CHAMELEON_Complex64_t * ) malloc( sizeof( CHAMELEON_Complex64_t ) * LDC * N  );
        C    = ( CHAMELEON_Complex64_t * ) malloc( sizeof( CHAMELEON_Complex64_t ) * LDC * N  );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descC,    C,    LDC );

    if ( rank == 0 ) {
        info_solution = check_zsymm_std( args, matrix_type, side, uplo,
                                         M, N, alpha, A, LDA, B, LDB,
                                         beta, Cref, C, LDC );
        free(A);
        free(B);
        free(C);
        free(Cref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed matrices rank k operations.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] N
 *          The order of the matrix C and number of rows of the matrices A and B.
 *
 * @param[in] K
 *          The number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] B
 *          The matrix B - only used for her2k and syr2k.
 *
 * @param[in] LDB
 *          The leading dimension of the matrix B.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] Cref
 *          The symmetric matrix Cref.
 *
 * @param[in] C
 *          The symmetric matrix of the computed result.
 *
 * @param[in] LDC
 *          The leading dimension of the matrices Cref and C.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsyrk_std( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans, int N, int K, CHAMELEON_Complex64_t alpha,
                     CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *B, int LDB, CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Cref,
                     CHAMELEON_Complex64_t *C, int LDC )
{
    int     info_solution = 0;
    double  Anorm, Bnorm, Crefnorm, Cchamnorm, Clapacknorm, Rnorm, result;
    double *work = malloc(sizeof(double)*N);

    Bnorm = 0.;
    if ( trans == ChamNoTrans ) {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, K, A, LDA );
        if ( B != NULL ) {
            Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', N, K, B, LDB );
        }
    }
    else {
        Anorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, N, A, LDA );
        if ( B != NULL ) {
            Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', K, N, B, LDB );
        }
    }

    /* Computes the norms for comparing */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Crefnorm  = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        Cchamnorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, C,    LDC, work );
    }
    else
#endif
    {
        Crefnorm  = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
        Cchamnorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, C,    LDC, work );
    }

    cham_fixdbl_t eps = testing_getaccuracy();
    double ABnorm;

    /* Makes the multiplication with the core function */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        if ( B == NULL ) {
            cblas_zherk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, creal(alpha), A, LDA, creal(beta), Cref, LDC );
            ABnorm = Anorm * Anorm;
        }
        else {
            cblas_zher2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, B, LDB, creal(beta), Cref, LDC );
            ABnorm = 2. * Anorm * Bnorm;
        }

        Clapacknorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
    }
    else
#endif
    {
        if ( B == NULL ) {
            cblas_zsyrk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), Cref, LDC );
            ABnorm = Anorm * Anorm;
        }
        else {
            cblas_zsyr2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                            N, K, CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC );
            ABnorm = 2. * Anorm * Bnorm;
        }

        Clapacknorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'I', chameleon_lapack_const(uplo), N, Cref, LDC, work );
    }

    CORE_ztradd( uplo, ChamNoTrans, N, N, -1., C, LDC, 1., Cref, LDC );

    /* Computes the norm with the core function's result */
#if defined(PRECISION_z) || defined(PRECISION_c)
    if ( matrix_type == ChamHermitian ) {
        Rnorm = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), N, Cref, LDC, NULL );
    }
    else
#endif
    {
        Rnorm = LAPACKE_zlansy_work( LAPACK_COL_MAJOR, 'M', chameleon_lapack_const(uplo), N, Cref, LDC, NULL );
    }
    result = Rnorm / ((ABnorm + Crefnorm) * K * eps);
    run_arg_add_double( args, "||A||", Anorm );
    run_arg_add_double( args, "||B||", Bnorm );
    run_arg_add_double( args, "||C||", Crefnorm );
    run_arg_add_double( args, "||R||", Rnorm );
    /* Verifies if the result is inside a threshold */
    if ( isinf(Clapacknorm) || isinf(Cchamnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(work);

    (void)matrix_type;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed matrix rank k operation with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] matrix_type
 *          Whether it is a general, triangular, hermitian or symmetric matrix.
 *
 * @param[in] uplo
 *          Whether it is a upper triangular matrix, a lower triangular matrix or a general matrix.
 *
 * @param[in] trans
 *          Whether the first product element is transposed, conjugate transposed or not transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descB
 *          The descriptor of the matrix B - only used for her2k and syr2k.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in] descCref
 *          The descriptor of the symmetric matrix C.
 *
 * @param[in] descC
 *          The descriptor of the symmetric matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zsyrk( run_arg_list_t *args, cham_mtxtype_t matrix_type, cham_uplo_t uplo, cham_trans_t trans,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descB,
                 CHAMELEON_Complex64_t beta, CHAM_desc_t *descCref, CHAM_desc_t *descC )
{
    int info_solution           = 0;
    int An, Bn, K, N;
    int LDA, LDB, LDC;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    CHAMELEON_Complex64_t *Cref = NULL;
    CHAMELEON_Complex64_t *C    = NULL;

    if ( trans == ChamNoTrans ) {
        N   = descA->m;
        K   = descA->n;
        An  = K;
        Bn  = K;
        LDA = N;
        LDB = N;
    }
    else {
        N   = descA->n;
        K   = descA->m;
        An  = N;
        Bn  = N;
        LDA = K;
        LDB = K;
    }
    LDC = N;

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDA * An  );
        if ( descB != NULL ) {
            B = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDB * Bn  );
        }
        Cref = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDC * N  );
        C    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDC * N  );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( uplo,           descCref, Cref, LDC );
    CHAMELEON_zDesc2Lap( uplo,           descC,    C,    LDC );
    if ( descB != NULL ) {
        CHAMELEON_zDesc2Lap( ChamUpperLower, descB, B, LDB );
    }

    if ( rank == 0 ) {

        info_solution = check_zsyrk_std( args, matrix_type, uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, Cref, C, LDC );

        free(A);
        free(C);
        free(Cref);
        if ( descB != NULL ) {
            free(B);
        }
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares two core function computed matrix triangular product.
 *
 *******************************************************************************
 *
 * @param[in] check_func
 *          Whether it is a triangular product or a triangular linear solution.
 *
 * @param[in] side
 *          Whether A appears on the left or on the right of the product.
 *
 * @param[in] uplo
 *          Whether A is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] trans
 *          Whether A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] diag
 *          Whether A is a unitary diagonal matrix or not.
 *
 * @param[in] M
 *          The number of rows of the matrix B and the order of the matrix A if side is left.
 *
 * @param[in] N
 *          The number of columns of the matrix B and the order of the matrix A if side is right.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A.
 *
 * @param[in] Bref
 *          The matrix Bref.
 *
 * @param[in] B
 *          The matrix of the computed result.
 *
 * @param[in] LDB
 *          The leading dimension of the matrices Bref and B.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrmm_std( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag, int M, int N,
                     CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA, CHAMELEON_Complex64_t *Bref, CHAMELEON_Complex64_t *B, int LDB )
{
    CHAMELEON_Complex64_t mzone = -1.0;
    int     info_solution, An;
    double  Anorm, Bnorm, Rnorm, result;
    char    normTypeA, normTypeB;
    double *work;
    cham_fixdbl_t eps = testing_getaccuracy();

    /* Computes the norms for comparing */
    if ( side == ChamLeft ) {
        normTypeA = 'O';
        if ( trans == ChamNoTrans ) {
            normTypeA = 'I';
        }
        normTypeB = 'O';
        An = M;
    }
    else {
        normTypeA = 'O';
        if ( trans != ChamNoTrans ) {
            normTypeA = 'I';
        }
        normTypeB = 'I';
        An = N;
    }

    work  = malloc( sizeof(double) * An );
    Anorm = LAPACKE_zlantr_work( LAPACK_COL_MAJOR, normTypeA,
                                 chameleon_lapack_const(uplo), chameleon_lapack_const(diag),
                                 An, An, A, LDA, work );
    free( work );

    Bnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, normTypeB, M, N, B, LDB );

    /* Makes the multiplication with the core function */
    if (check_func == CHECK_TRMM) {
        cblas_ztrmm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                    (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB);
    }
    else {
        cblas_ztrsm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                    (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB);
    }

    /* Computes the norm with the core function's result */
    cblas_zaxpy( LDB * N, CBLAS_SADDR(mzone), Bref, 1, B, 1 );
    Rnorm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, B, LDB, NULL );

    result = Rnorm / ((Anorm + Bnorm) * An * eps);

    /* Verifies if the result is inside a threshold */
    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    (void)Bnorm;
    (void)args;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Compares a Chameleon computed matrix triangular product with a core function computed one.
 *
 *******************************************************************************
 *
 * @param[in] check_func
 *          Whether it is a triangular product or a triangular linear solution.
 *
 * @param[in] side
 *          Whether A appears on the left or on the right of the product.
 *
 * @param[in] uplo
 *          Whether A is a upper triangular matrix or a lower triangular matrix.
 *
 * @param[in] trans
 *          Whether A is transposed, conjugate transposed or not transposed.
 *
 * @param[in] diag
 *          Whether A is a unitary diagonal matrix or not.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] descA
 *          The descriptor of the matrix A.
 *
 * @param[in] descBref
 *          The descriptor of the matrix Bref.
 *
 * @param[in] descB
 *          The descriptor of the matrix of the Chameleon computed result.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_ztrmm( run_arg_list_t *args, int check_func, cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
                 CHAMELEON_Complex64_t alpha, CHAM_desc_t *descA, CHAM_desc_t *descBref, CHAM_desc_t *descB )
{
    int info_solution           = 0;
    int M                       = descBref->m;
    int N                       = descBref->n;
    int rank                    = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *A    = NULL;
    CHAMELEON_Complex64_t *Bref = NULL;
    CHAMELEON_Complex64_t *B    = NULL;
    int An                      = ( side == ChamLeft )? M : N;
    int LDA                     = An;
    int LDB                     = M;

    if ( rank == 0 ) {
        A    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDA * An );
        Bref = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDB * N  );
        B    = (CHAMELEON_Complex64_t *)malloc( sizeof(CHAMELEON_Complex64_t) * LDB * N  );
    }

    /* Creates the LAPACK version of the matrices */
    CHAMELEON_zDesc2Lap( uplo,           descA,    A,    LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descB,    B,    LDB );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descBref, Bref, LDB );

    if ( rank == 0 ) {

        info_solution = check_ztrmm_std( args, check_func, side, uplo, trans, diag, M, N, alpha, A, LDA, Bref, B, LDB );

        free(A);
        free(B);
        free(Bref);
    }

    /* Broadcasts the result from the main processus */
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return info_solution;
}

#endif /* defined(CHAMELEON_SIMULATION) */
