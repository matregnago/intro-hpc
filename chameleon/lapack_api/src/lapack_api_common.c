/**
 *
 * @file lapack_api_common.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api common functions
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-11-14
 *
 */
#include "lapack_api_common.h"

/**
 * @brief Convert the input char BLAS trans parameter to a compatible parameter
 * for the Cblas API.
 * @param[in] trans The input char BLAS trans parameter
 * @return The CBLAS equivalent parameter (CblasNoTrans, CblasTrans or
 * CblasConjTrans).
 */
int chameleon_blastocblas_trans(const char *trans)
{
    if ( (*trans == 'N') || (*trans == 'n') ) {
        return CblasNoTrans;
    } else if ( (*trans == 'T') || (*trans == 't') ) {
        return CblasTrans;
    } else if ( (*trans == 'C') || (*trans == 'c') ) {
        return CblasConjTrans;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_blastocblas_trans", "illegal value of BLAS transpose parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}

/**
 * @brief Convert the input char BLAS side parameter to a compatible parameter
 * for the Cblas API.
 * @param[in] uplo The input char BLAS side parameter
 * @return The CBLAS equivalent parameter (CblasLeft or CblasRight).
 */
int chameleon_blastocblas_side(const char *side)
{
    if ( (*side == 'L') || (*side == 'l') ) {
        return CblasLeft;
    } else if ( (*side == 'R') || (*side == 'r') ) {
        return CblasRight;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_blastocblas_side", "illegal value of BLAS side parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}

/**
 * @brief Convert the input char BLAS uplo parameter to a compatible parameter
 * for the Cblas API.
 * @param[in] uplo The input char BLAS uplo parameter
 * @return The CBLAS equivalent parameter (CblasUpper or CblasLower).
 */
int chameleon_blastocblas_uplo(const char *uplo)
{
    if ( (*uplo == 'U') || (*uplo == 'u') ) {
        return CblasUpper;
    } else if ( (*uplo == 'L') || (*uplo == 'l') ) {
        return CblasLower;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_blastocblas_uplo", "illegal value of BLAS uplo parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}

/**
 * @brief Convert the input char BLAS diag parameter to a compatible parameter
 * for the Cblas API.
 * @param[in] diag The input char BLAS diag parameter
 * @return The CBLAS equivalent parameter (CblasUnit or CblasNonUnit).
 */
int chameleon_blastocblas_diag(const char *diag)
{
    if ( (*diag == 'U') || (*diag == 'u') ) {
        return CblasUnit;
    } else if ( (*diag == 'N') || (*diag == 'n') ) {
        return CblasNonUnit;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_blastocblas_diag", "illegal value of BLAS diag parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}

/**
 * @brief Convert the input char LAPACK norm parameter to a compatible parameter
 * for the Chameleon API.
 * @param[in] norm The input char LAPACK norm parameter
 * @return The Chameleon equivalent parameter (ChamMaxNorm, ChamOneNorm, etc).
 */
int chameleon_lapacktochameleon_norm(const char *norm)
{
    if ( (*norm == 'M') || (*norm == 'm') ) {
        return ChamMaxNorm;
    } else if ( (*norm == '1') || (*norm == 'O') || (*norm == 'o') ) {
        return ChamOneNorm;
    } else if ( (*norm == 'I') || (*norm == 'i') ) {
        return ChamInfNorm;
    } else if ( (*norm == 'F') || (*norm == 'f') || (*norm == 'E') || (*norm == 'e') ) {
        return ChamFrobeniusNorm;
    } else {
        fprintf(stderr, "CHAMELEON ERROR: %s(): %s\n", "chameleon_lapacktochameleon_norm", "illegal value of LAPACK norm parameter");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}