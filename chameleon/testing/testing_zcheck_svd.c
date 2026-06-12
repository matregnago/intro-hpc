/**
 *
 * @file testing_zcheck_svd.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @date 2025-12-19
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
 * @brief Checks if the Chameleon SVD algorithm works: Ainit = U * mat( S ) * Vt.
 *        - U and Vt should be orthogonal.
 *        - Sinit and S should be the same.
 *        - Ainit = U * mat( S ) * Vt.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V^H.
 *
 * @param[in] M
 *          The number of rows of the matrices Ainit, A and U.
 *
 * @param[in] N
 *          The number of columns of the matrices Ainit, A and Vt.
 *
 * @param[in] Ainit
 *          The matrix Ainit (initial matrix A).
 *
 * @param[in] A
 *          The matrix A after the SVD, can contain parts of the matrix U or Vt
 *          or nothing (depending on the values of jobu and jobvt).
 *
 * @param[in] LDA
 *          The leading dimension of the matrices A and Ainit.
 *
 * @param[in] Sinit
 *          The vector with the singular values of the matrix Ainit
 *          (contains the K = min(M, N) singular values of Ainit).
 *
 * @param[in] S
 *          The vector with the singular values of the matrix Ainit
 *          computed by the Chameleon SVD algorithm.
 *
 * @param[in] U
 *          The orthogonal matrix U computed by the Chameleon SVD algorithm can
 *          contain all of U, a part of U or nothing (NULL) depending on the value of jobu;
 *          if jobu == AllVec : dimension: M * M
 *          if jobu == SVec :   dimension: M * min(M, N)
 *          if jobu == NoVec or OVec : U = NULL
 *
 * @param[in] LDU
 *          The leading dimension of the matrix U.
 *
 * @param[in] Vt
 *          The orthogonal matrix Vt computed by the Chameleon SVD algorithm can
 *          contain all of Vt, a part of Vt or nothing (NULL) depending on the value of jobvt;
 *          if jobuvt == AllVec : dimension: N * N
 *          if jobvt == SVec :    dimension: min(M, N) * N
 *          if jobvt == NoVec or OVec : Vt = NULL
 *
 * @param[in] LDVt
 *          The leading dimension of the matrix Vt.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgesvd_std( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, int M, int N, CHAMELEON_Complex64_t *Ainit, CHAMELEON_Complex64_t *A, int LDA,
                     double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt )
{
    int info_solution = 0;
    double result;
    int k;
    int K = chameleon_min(M, N);
    cham_fixdbl_t eps = testing_getaccuracy();

    /* Checks if U is orthogonal */
    switch ( jobu ) {
    case ChamAllVec:
        info_solution += check_zortho_std( args, M, M, U, LDU );
        break;
    case ChamSVec:
        info_solution += check_zortho_std( args, M, K, U, LDU );
        break;
    case ChamOVec:
        info_solution += check_zortho_std( args, M, K, A, LDA );
        break;
    default:
       ;
    }

    /* Checks if Vt is orthogonal */
    switch ( jobvt ) {
    case ChamAllVec:
        info_solution += check_zortho_std( args, N, N, Vt, LDVt );
        break;
    case ChamSVec:
        info_solution += check_zortho_std( args, K, N, Vt, LDVt );
        break;
    case ChamOVec:
        info_solution += check_zortho_std( args, K, N, A, LDA );
        break;
    default:
       ;
    }

    /* checks if Sinit and S are the same */
    double maxSV     = chameleon_max( fabs(Sinit[0]), fabs(S[0]) );
    double maxSVdiff = fabs( fabs(Sinit[0]) - fabs(S[0]) );
    double maxtmp, maxdifftmp;

    for ( k = 1; k < K; k++ ) {
        maxdifftmp = fabs( fabs(Sinit[k]) - fabs(S[k]) );
        maxtmp     = chameleon_max( fabs(Sinit[k]), fabs(S[k]) );

        /* Update */
        maxSV     = chameleon_max( maxtmp, maxSV );
        maxSVdiff = chameleon_max( maxdifftmp, maxSVdiff );
    }

    result = maxSVdiff / ( K * eps * maxSV );

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        info_solution += 1;
    }

    if ( (jobu == ChamAllVec) && (jobvt == ChamAllVec) ) {
        /* To do: create mat( S ) */
    }
    result = info_solution;

    run_arg_add_double( args, "||R||", result );

    (void)Ainit;
    return info_solution;
}

/**
 ********************************************************************************
 *
 * @ingroup testing
 *
 * @brief Checks if the Chameleon SVD algorithm works: descAinit = U * mat( S ) * Vt.
 *        - U and Vt should be orthogonal.
 *        - Sinit and S should be the same.
 *        - descAinit = U * mat( S ) * Vt.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V^H.
 *
 * @param[in] descAinit
 *          The descriptor of the matrix Ainit (initial matrix A).
 *
 * @param[in] descA
 *          The descriptor of the matrix A after the SVD, can contain parts of the matrix
 *          U or Vt or nothing (depending on the values of jobu and jobvt).
 *
 * @param[in] Sinit
 *          The vector with the singular values of the matrix Ainit
 *          (contains the K = min(M, N) singular values of Ainit).
 *
 * @param[in] S
 *          The vector with the singular values of the matrix Ainit
 *          computed by the Chameleon SVD algorithm.
 *
 * @param[in] U
 *          The orthogonal matrix U computed by the Chameleon SVD algorithm can
 *          contain all of U, a part of U or nothing (NULL) depending on the value of jobu;
 *          if jobu == AllVec : dimension: M * M
 *          if jobu == SVec :   dimension: M * min(M, N)
 *          if jobu == NoVec or OVec : U = NULL
 *
 * @param[in] LDU
 *          The leading dimension of the matrix U.
 *
 * @param[in] Vt
 *          The orthogonal matrix Vt computed by the Chameleon SVD algorithm can
 *          contain all of Vt, a part of Vt or nothing (NULL) depending on the value of jobvt;
 *          if jobuvt == AllVec : dimension: N * N
 *          if jobvt == SVec :    dimension: min(M, N) * N
 *          if jobvt == NoVec or OVec : Vt = NULL
 *
 * @param[in] LDVt
 *          The leading dimension of the matrix Vt.
 *
 * @retval 0 successfull comparison
 *
 *******************************************************************************
 */
int check_zgesvd( run_arg_list_t *args, cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *descAinit, CHAM_desc_t *descA,
                  double *Sinit, double *S, CHAMELEON_Complex64_t *U, int LDU, CHAMELEON_Complex64_t *Vt, int LDVt )
{
    int info_solution = 0;
    int rank = CHAMELEON_Comm_rank();
    CHAMELEON_Complex64_t *Ainit = NULL;
    CHAMELEON_Complex64_t *A     = NULL;
    int M    = descA->m;
    int N    = descA->n;
    int LDA  = descA->lm;

    if ( rank == 0 ) {
        Ainit = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        A     = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    }

    CHAMELEON_zDesc2Lap( ChamUpperLower, descAinit, Ainit, LDA );
    CHAMELEON_zDesc2Lap( ChamUpperLower, descA,     A,     LDA );

    if ( rank == 0 ) {

        info_solution = check_zgesvd_std( args, jobu, jobvt, M, N, Ainit, A, LDA, Sinit, S, U, LDU, Vt, LDVt  );

        free( Ainit );
        free( A     );
    }

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

    return info_solution;
}

#endif /* defined(CHAMELEON_SIMULATION) */
