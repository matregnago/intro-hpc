/**
 *
 * @file hmat_z.h
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2019-2019 Universidad Jaume I. All rights reserved.
 *
 * @brief Chameleon CPU kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.4.0
 * @author Rocio Carratala-Saez
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#ifndef _hmat_z_h_
#define _hmat_z_h_

HMAT_API int hmat_zgetrf( hmat_matrix_t *A );
HMAT_API int hmat_zpotrf( hmat_matrix_t *A );
HMAT_API int hmat_zgemm( char transA, char transB, void* alpha, hmat_matrix_t* A, hmat_matrix_t* B, void* beta, hmat_matrix_t* C );
HMAT_API int hmat_ztrsm( char side, char uplo, char trans, char diag, int m, int n,
			 void* alpha, hmat_matrix_t* A, int is_b_hmat, void* B );
HMAT_API int hmat_zgemv( char transA, void* alpha, hmat_matrix_t* A,
			 void* B, void* beta, void* C, int nrhs );
HMAT_API hmat_matrix_t *hmat_zread( void *buffer );
HMAT_API size_t         hmat_zsize( hmat_matrix_t *hmat );
HMAT_API void           hmat_zwrite( hmat_matrix_t *hmat, char *ptr );

HMAT_API void hmat_zdestroy( hmat_matrix_t *hmat );

#endif /* _hmat_z_h_ */
