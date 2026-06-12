/**
 *
 * @file hmat_z.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2019-2019 Universidad Jaume I. All rights reserved.
 *
 * @brief Chameleon interface for H-Mat kernels
 *
 * @version 1.4.0
 * @author Rocio Carratala-Saez
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/hmat.h"

/**
 * @brief The hmat interface to the C++ functions
 */
hmat_interface_t hmat_zinterface;

/**
 * @brief Constructor of the C++ interface
 */
void __hmat_zinit() __attribute__(( constructor ));
void __hmat_zinit() {
    hmat_init_default_interface( &hmat_zinterface, HMAT_DOUBLE_COMPLEX );
}

/**
 * @brief Destructor of the C++ interface
 */
void __hmat_zfini() __attribute__(( destructor ));
void __hmat_zfini() {
    hmat_zinterface.finalize();
}

/**
 * @brief Cholesky factorization of an H-Matrix
 *
 * @param[in,out] A
 *      On entry, the input matrix A.
 *      On exit, the factorized matrix A.
 *
 * @return 0 (Return an integer to match prototype of LAPACK)
 */
int hmat_zpotrf( hmat_matrix_t *A ) {
    hmat_factorization_context_t ctx_facto;
    hmat_factorization_context_init( &ctx_facto );
    ctx_facto.factorization = hmat_factorization_llt;
    hmat_zinterface.factorize_generic( A, &ctx_facto );
    return 0;
}

/**
 * @brief LU factorization of an H-Matrix
 *
 * @param[in,out] A
 *      On entry, the input matrix A.
 *      On exit, the factorized matrix A.
 *
 * @return 0 (Return an integer to match prototype of LAPACK)
 */
int hmat_zgetrf( hmat_matrix_t *A ) {
    hmat_factorization_context_t ctx_facto;
    hmat_factorization_context_init( &ctx_facto );
    ctx_facto.factorization = hmat_factorization_lu;
    hmat_zinterface.factorize_generic( A, &ctx_facto );
    return 0;
}

/**
 * @brief Matrix-Matrix product of H-Matrices
 *
 * Computes \f[ C = \alpha op(A) op(B) + \beta C \f]
 *
 * where op(A) = A, A^t or conj(A^t)
 *
 * @param[in] transA
 *      'N' : op(A) = A
 *      'T' : op(A) = A^t
 *      'C' : op(A) = conj(A^t)
 *
 * @param[in] transB
 *      'N' : op(B) = B
 *      'T' : op(B) = B'
 *      'C' : op(B) = conj(B')
 *
 * @param[in] alpha
 *      The scalar alpha.
 *
 * @param[in] alpha
 *      The scalar alpha.
 *
 * @param[in] A
 *      On entry, the input matrix A.
 *
 * @param[in] B
 *      On entry, the input matrix B.
 *
 * @param[in] beta
 *      The scalar beta.
 *
 * @param[in,out] C
 *      On entry, the input matrix C.
 *      On exit, the updated matrix C.
 *
 * @return 0 (Return an integer to match prototype of LAPACK)
 */
int hmat_zgemm( char transA, char transB, void* alpha, hmat_matrix_t* A,
		hmat_matrix_t* B, void* beta, hmat_matrix_t* C ) {
    return hmat_zinterface.gemm( transA, transB, alpha, A, B, beta, C );
}

/**
 * @brief Matrix-Matrix product of an H-Matrix with a full-rank matrix.
 *
 * Computes \f[ C = \alpha op(A) B + \beta C \f]
 *
 * where op(A) = A, A^t or conj(A^t), and A is an H-Matrix, and B and C are
 * full-rank matrices.
 *
 * @param[in] transA
 *      'N' : op(A) = A
 *      'T' : op(A) = A^t
 *      'C' : op(A) = conj(A^t)
 *
 * @param[in] alpha
 *      The scalar alpha.
 *
 * @param[in] alpha
 *      The scalar alpha.
 *
 * @param[in] A
 *      On entry, the input H-Matrix A.
 *
 * @param[in] B
 *      On entry, the full-rank matrix B.
 *
 * @param[in] beta
 *      The scalar beta.
 *
 * @param[in,out] C
 *      On entry, the full-rank matrix C.
 *      On exit, the updated full-rank matrix C.
 *
 * @return 0 on success, -1 otherwise
 */
int hmat_zgemv( char transA, void* alpha, hmat_matrix_t* A,
		void* B, void* beta, void* C, int nrhs ) {
    return hmat_zinterface.gemm_scalar( transA, alpha, A, B, beta, C, nrhs );
}

/**
 * @brief Solve a triangular system \f[ A x = b \f] with A being an H-Matrix
 *
 * Solves \f[ \alpha op(A) X = B \f] or  \f[ \alpha X op(A) = B \f]
 *
 * where op(A) = A, A^t or conj(A^t), and A is an H-Matrix, and B and C are
 * full-rank matrices.
 *
 * @param[in] side
 *      'L' : Solves \alpha A X = B
 *      'R' : Solves \alpha X A = B
 *
 * @param[in] uplo
 *      'L' : A is lower triangular, upper part is not referenced
 *      'U' : A is upper triangular, lower part is not referenced
 *
 * @param[in] trans
 *      'N' : op(A) = A
 *      'T' : op(A) = A^t
 *      'C' : op(A) = conj(A^t)
 *
 * @param[in] uplo
 *      'U' : A has a unitary diagonal which is ot referenced
 *      'N' : A has a non unitary diagonal
 *
 * @param[in] m
 *      The number of rows of B
 *
 * @param[in] n
 *      The number of columns of B
 *
 * @param[in] alpha
 *      The scalar alpha.
 *
 * @param[in] A
 *      On entry, the input H-Matrix A.
 *
 * @param[in] is_b_hmat
 *      if true, B is an H-Matrix, otherwise B is full-rank and stored in lapack layout
 *
 * @param[in] B
 *      On entry, the H-Mat of full-rank matrix B.
 *      On entry, the solution of the trinagular system.
 *
 * @return 0 on success, -1 otherwise
 */
int hmat_ztrsm( char side, char uplo, char trans, char diag, int m, int n,
		void* alpha, hmat_matrix_t* A, int is_b_hmat, void* B ) {
    return hmat_zinterface.trsm( side, uplo, trans, diag, m, n, alpha, A, is_b_hmat, B );
}

/**
 * @brief Read/Unpack an H-Matrix from a packed format into a contiguous buffer
 *
 * @param[in] buffer
 *      The buffer that contains the packed H-Matrix.
 *
 * @return The unpack H-Matrix
 */
hmat_matrix_t *hmat_zread( void *buffer ) {
    hmat_buffer_comm_t buffer_struct = { 0, (char*)buffer };
    hmat_matrix_t *hmat = hmat_zinterface.read_struct( (hmat_iostream)(&buffer_comm_read),
                                                       (void *)(&buffer_struct) );
    hmat_zinterface.read_data( hmat, (hmat_iostream)(&buffer_comm_read), (void *)(&buffer_struct) );
    return hmat;
}

/**
 * @brief Compute the size required to store the H-Matrix in a packed format
 *
 * @param[in] hmat
 *      The H-Matrix for which the size must be computed
 *
 * @return The size of the H-Matrix (structure+data)
 */
size_t hmat_zsize( hmat_matrix_t *hmat ) {
    size_t size = 0;
    hmat_zinterface.write_struct( hmat, (hmat_iostream)(&buffer_comm_size), (void *)(&size) );
    hmat_zinterface.write_data(   hmat, (hmat_iostream)(&buffer_comm_size), (void *)(&size) );
    return size;
}

/**
 * @brief Write/pack an H-Matrix in a packed format in a given buffer
 *
 * @param[in] hmat
 *      The H-Matrix to pack
 *
 * @param[in,out] ptr
 *      On entry the allocated buffer with sufficent space to store the H-Matrix.
 *      On exit, contains the packed H-Matrix.
 */
void hmat_zwrite( hmat_matrix_t *hmat, char *ptr ) {
    hmat_buffer_comm_t buffer_struct = { 0, ptr };
    hmat_zinterface.write_struct( hmat, (hmat_iostream)(&buffer_comm_write), (void *)(&buffer_struct) );
    hmat_zinterface.write_data(   hmat, (hmat_iostream)(&buffer_comm_write), (void *)(&buffer_struct) );
}

/**
 * @brief Free the data structure (structure+data) associated to an H-Matrix
 *
 * @param[in] hmat
 *      The H-Matrix to destroy
 */
void hmat_zdestroy( hmat_matrix_t *hmat ) {
    hmat_zinterface.destroy( hmat );
}
