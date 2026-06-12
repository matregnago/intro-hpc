/**
 *
 * @file gpucublas_ztile.h
 *
 * @copyright 2025-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CUDA tile-based kernel interface header
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @author Florent Pruvost
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpucublas_ztile_h_
#define _gpucublas_ztile_h_

#include "chameleon/struct.h"
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuComplex.h>

/**
 *  Declarations of CUDA tile-based kernels - alphabetical order
 */
int TCUDA_zgeadd( cham_trans_t trans, int m, int n, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const cuDoubleComplex *beta, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_zgemerge( cham_side_t side, cham_diag_t diag, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_zgemm( cham_trans_t transA, cham_trans_t transB, int m, int n, int k, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, const cuDoubleComplex *beta, CHAM_tile_t *C, cublasHandle_t handle );
#if defined( PRECISION_z ) || defined( PRECISION_c )
int TCUDA_zhemm( cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, const cuDoubleComplex *beta, CHAM_tile_t *C, cublasHandle_t handle );
int TCUDA_zher2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, const double *beta, CHAM_tile_t *C, cublasHandle_t handle );
int TCUDA_zherk( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const CHAM_tile_t *A, const double *beta, CHAM_tile_t *C, cublasHandle_t handle );
#endif
int TCUDA_zherfb( cham_uplo_t uplo, int n, int k, int ib, int nb, const CHAM_tile_t *A, const CHAM_tile_t *T, CHAM_tile_t *C, CHAM_tile_t *WORK, int ldwork, cublasHandle_t handle );
int TCUDA_zlarfb( cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *C, CHAM_tile_t *WORK, int ldwork, cublasHandle_t handle );
int TCUDA_zlacpy( cham_uplo_t uplo, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_zlacpyx( cham_uplo_t uplo, int M, int N, int displA, const CHAM_tile_t *A, int LDA, int displB, CHAM_tile_t *B, int LDB, cublasHandle_t handle );
int TCUDA_zlaset( cham_uplo_t uplo, int n1, int n2, const cuDoubleComplex *alpha, const cuDoubleComplex *beta, CHAM_tile_t *A, cublasHandle_t handle );
int TCUDA_zlatro( cham_uplo_t uplo, cham_trans_t trans, int M, int N, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_zparfb( cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_zpotrf( cham_uplo_t uplo, int n, CHAM_tile_t *A, cuDoubleComplex *WORK, int lwork, int *d_info, cusolverDnHandle_t handle );
int TCUDA_zsymm( cham_side_t side, cham_uplo_t uplo, int m, int n, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, const cuDoubleComplex *beta, CHAM_tile_t *C, cublasHandle_t handle );
int TCUDA_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const CHAM_tile_t *B, const cuDoubleComplex *beta, CHAM_tile_t *C, cublasHandle_t handle );
int TCUDA_zsyrk( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const cuDoubleComplex *alpha, const CHAM_tile_t *A, const cuDoubleComplex *beta, CHAM_tile_t *C, cublasHandle_t handle );
int TCUDA_ztpmqrt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *A, CHAM_tile_t *B, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_ztpmlqt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *A, CHAM_tile_t *B, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_ztrmm( cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int m, int n, const cuDoubleComplex *alpha, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_ztrsm( cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag, int m, int n, const cuDoubleComplex *alpha, const CHAM_tile_t *A, CHAM_tile_t *B, cublasHandle_t handle );
int TCUDA_ztsmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_zttmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_ztsmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_zttmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, CHAM_tile_t *A1, CHAM_tile_t *A2, const CHAM_tile_t *V, const CHAM_tile_t *T, CHAM_tile_t *WORK, int lwork, cublasHandle_t handle );
int TCUDA_zunmlqt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAM_tile_t *A, const CHAM_tile_t *T, CHAM_tile_t *C, CHAM_tile_t *WORK, int ldwork, cublasHandle_t handle );
int TCUDA_zunmqrt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const CHAM_tile_t *A, const CHAM_tile_t *T, CHAM_tile_t *C, CHAM_tile_t *WORK, int ldwork, cublasHandle_t handle );

#endif /* _gpucublas_ztile_h_ */
