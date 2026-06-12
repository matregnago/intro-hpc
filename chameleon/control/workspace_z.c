/**
 *
 * @file workspace_z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon precision dependent workspace routines
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "control/common.h"
#include "control/workspace.h"

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgeev - Allocates workspace for CHAMELEON_zgeev or
 *  CHAMELEON_zgeev_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile Hessenberg.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgeev(int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(N, N, CHAMELEON_FUNC_ZGEEV, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgehrd - Allocates workspace for CHAMELEON_zgehrd or
 *  CHAMELEON_zgehrd_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile Hessenberg.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgehrd(int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(N, N, CHAMELEON_FUNC_ZGEHRD, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgebrd - Allocates workspace for CHAMELEON_zgebrd or CHAMELEON_zgebrd_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile BRD.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgebrd(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGEBRD, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgels - Allocates workspace for CHAMELEON_zgels or
 *  CHAMELEON_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile QR or the tile LQ factorization.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgels(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgels_Tile - Allocates tile workspace for
 *  CHAMELEON_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, CHAM_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgels_Tile(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgeqrf - Allocates workspace for CHAMELEON_zgeqrf or
 *  CHAMELEON_zgeqrf_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors
 *          required by the tile QR factorization.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgeqrf(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgeqrf_Tile - Allocates tile workspace for
 *  CHAMELEON_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, CHAM_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgeqrf_Tile(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgelqf - Allocates workspace for CHAMELEON_zgelqf or
 *  CHAMELEON_zgelqf_Tile routines.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile LQ
 *          factorization.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgelqf(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgelqf_Tile - Allocates tile workspace for CHAMELEON_zgels_Tile routine.
 *
 ******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, CHAM_desc_t *on workspace handle for storage of the extra
 *          T factors required by the tile QR or the tile LQ factorization.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgelqf_Tile(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGELS, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgesv - Allocates workspace for CHAMELEON_zgesv or
 *  CHAMELEON_zgesv_Tile routines.
 *
 ******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A.
 *          N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors
 *          required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required
 *          by the tile LU factorization (not equivalent to LAPACK).
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgesv_incpiv(int N, CHAM_desc_t **descL, int **IPIV, int p, int q) {
    return chameleon_alloc_ipiv(N, N, CHAMELEON_FUNC_ZGESV, ChamComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgesv_Tile - Allocates workspace for CHAMELEON_zgesv_Tile
 *  routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A.
 *          N >= 0.
 *
 * @param[out] descL
 *          On exit, CHAMELEON descriptor on workspace handle for storage of the
 *          extra L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by
 *          the tile LU factorization (not equivalent to LAPACK).
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile(int N, CHAM_desc_t **descL, int **IPIV, int p, int q)
{
    return chameleon_alloc_ipiv(N, N, CHAMELEON_FUNC_ZGESV, ChamComplexDouble, descL, (void**)IPIV, p, q);
}
/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgesvd - Allocates workspace for CHAMELEON_zgesvd or
 *  CHAMELEON_zgesvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgesvd(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZGESVD, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgetrf_incpiv - Allocates workspace for
 *  CHAMELEON_zgetrf_incpiv or CHAMELEON_zgetrf_incpiv_Tile or
 *  CHAMELEON_zgetrf_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 ******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_incpiv
 * @sa CHAMELEON_zgetrf_incpiv_Tile
 * @sa CHAMELEON_zgetrf_incpiv_Tile_Async
 *
 */
int CHAMELEON_Alloc_Workspace_zgetrf_incpiv(int M, int N, CHAM_desc_t **descL, int **IPIV, int p, int q) {
    return chameleon_alloc_ipiv(M, N, CHAMELEON_FUNC_ZGESV, ChamComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile - Allocates workspace for
 *  CHAMELEON_zgesv_incpiv_Tile or CHAMELEON_zgesv_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, CHAMELEON descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile(int N, CHAM_desc_t **descL, int **IPIV, int p, int q) {
    return chameleon_alloc_ipiv(N, N, CHAMELEON_FUNC_ZGESV, ChamComplexDouble, descL, (void**)IPIV, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zheev - Allocates workspace for CHAMELEON_zheev or CHAMELEON_zheev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zheev(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZHEEV, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zheevd - Allocates workspace for CHAMELEON_zheevd or CHAMELEON_zheevd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zheevd(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZHEEVD, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zhegv - Allocates workspace for CHAMELEON_zhegv or CHAMELEON_zhegv_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zhegv(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZHEGV, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zhegvd - Allocates workspace for CHAMELEON_zhegvd or CHAMELEON_zhegvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zhegvd(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZHEGVD, ChamComplexDouble, descT, p, q); }

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Alloc_Workspace_zhetrd - Allocates workspace for CHAMELEON_zhetrd or CHAMELEON_zhetrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 * @param[in] p
 *          Number of processes rows for the 2D block-cyclic distribution.
 *
 * @param[in] q
 *          Number of processes columns for the 2D block-cyclic distribution.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Alloc_Workspace_zhetrd(int M, int N, CHAM_desc_t **descT, int p, int q) {
    return chameleon_alloc_ibnb_tile(M, N, CHAMELEON_FUNC_ZHETRD, ChamComplexDouble, descT, p, q); }
