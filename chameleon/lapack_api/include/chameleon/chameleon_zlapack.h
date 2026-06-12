/**
 *
 * @file chameleon_zlapack.h
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api functions
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-01-22
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_zlapack_h_
#define _chameleon_zlapack_h_

#include <chameleon.h>
#include <coreblas/cblas_wrapper.h>

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout, Cblas/Lapacke interface) - alphabetical order
 */
void CHAMELEON_cblas_zgemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                            int M, int N, int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                               const CHAMELEON_Complex64_t *B, int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zhemm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                               const CHAMELEON_Complex64_t *B, int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zher2k( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                             int N, int K,
                             const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                const CHAMELEON_Complex64_t *B, int ldb,
                             double beta,             CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zherk( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                            int N, int K,
                            double alpha, const CHAMELEON_Complex64_t *A, int lda,
                            double beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zsymm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                               const CHAMELEON_Complex64_t *B, int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zsyr2k( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                             int N, int K,
                             const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                const CHAMELEON_Complex64_t *B, int ldb,
                             const void *beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_zsyrk( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                            int N, int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                            const void *beta,        CHAMELEON_Complex64_t *C, int ldc );

void CHAMELEON_cblas_ztrmm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            CBLAS_TRANSPOSE trans, CBLAS_DIAG diag,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                     CHAMELEON_Complex64_t *B, int ldb );

void CHAMELEON_cblas_ztrsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo,
                            CBLAS_TRANSPOSE trans, CBLAS_DIAG diag,
                            int M, int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, int lda,
                                                     CHAMELEON_Complex64_t *B, int ldb );

int CHAMELEON_lapacke_zlacpy( int matrix_layout, char uplo, int M, int N,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb );

double CHAMELEON_lapacke_zlange( int matrix_layout, char norm, int M, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlanhe( int matrix_layout, char norm, char uplo, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlansy( int matrix_layout, char norm, char uplo, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlantr( int matrix_layout, char norm, char uplo, char diag,
                                 int M, int N, const CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zlaset( int matrix_layout, char uplo, int M, int N,
                              CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zlauum( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zposv( int matrix_layout, char uplo, int N, int NRHS,
                             CHAMELEON_Complex64_t *A, int lda,
                             CHAMELEON_Complex64_t *B, int ldb );

int CHAMELEON_lapacke_zpotrf( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zpotri( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zpotrs( int matrix_layout, char uplo, int N, int NRHS,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb );

int CHAMELEON_lapacke_ztrtri( int matrix_layout, char uplo, char diag, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zgetrf( int matrix_layout, int M, int N,
                              CHAMELEON_Complex64_t *A, int lda, int *IPIV );

END_C_DECLS

#endif /* _chameleon_zlapack_h_ */
