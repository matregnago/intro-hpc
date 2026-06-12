/**
 *
 * @file starpu/codelet_zgetrf_nopiv_percol.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon getrf codelets to factorize the panel with no pivoting
 *
 * @version 1.4.0
 *
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

CHAMELEON_CL_CB( zgetrf_nopiv_percol_diag, cti_handle_get_m(task->handles[0]), 0, 0, M )
CHAMELEON_CL_CB( zgetrf_nopiv_percol_trsm, cti_handle_get_m(task->handles[0]), 0, 0, M )

#if !defined(CHAMELEON_SIMULATION)
static const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
static const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

static void cl_zgetrf_nopiv_percol_diag_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *tileA, *tileU;
    int                    m, n, k, lda, iinfo;
    RUNTIME_sequence_t    *sequence;
    RUNTIME_request_t     *request;
    CHAMELEON_Complex64_t *A, *row, pivot;

    tileA = cti_interface_get( descr[0] );
    tileU = cti_interface_get( descr[1] );

    starpu_codelet_unpack_args( cl_arg, &m, &n, &k, &iinfo, &sequence, &request );

    A   = tileA->mat;
    lda = tileA->ld;
    row = tileU->mat;

    /* Shift to the diagonal element */
    A += k * (lda + 1);

    /* Extract row into buffer */
    cblas_zcopy( n-k, A, lda, row, 1 );

    /* Perform update on current diagonal block directly here */
    if ( *row == 0. ) {
        if ( sequence->status == CHAMELEON_SUCCESS ) {
            RUNTIME_sequence_flush( NULL, sequence, request, iinfo+k+1 );
        }
        return;
    }

    pivot = 1. / *row;
    cblas_zscal( m-k-1, CBLAS_SADDR( pivot ), A + 1, 1 );

    CORE_zgemm( ChamNoTrans, ChamNoTrans,
                m-k-1, n-k-1, 1,
                mzone, A   + 1,       lda,
                       row + 1,       1,
                zone,  A   + 1 + lda, lda );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_nopiv_percol_diag, cl_zgetrf_nopiv_percol_diag_cpu_func )

void INSERT_TASK_zgetrf_nopiv_percol_diag( const RUNTIME_option_t *options,
                                                 int m, int n, int k,
                                                 const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *U, int Um, int Un,
                                                 int iinfo )
{
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_ACCESS_W( U, Um, Un );
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgetrf_nopiv_percol_diag,
        STARPU_VALUE, &m,                   sizeof(int),
        STARPU_VALUE, &n,                   sizeof(int),
        STARPU_VALUE, &k,                   sizeof(int),
        STARPU_VALUE, &iinfo,               sizeof(int),
        STARPU_VALUE, &(options->sequence), sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE, &(options->request),  sizeof(RUNTIME_request_t*),

        /* Task handles */
        STARPU_RW, RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,  RTBLKADDR(U, ChamComplexDouble, Um, Un),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgetrf_nopiv_percol_diag ),
        0 );
}

/*
 * Update column blocs
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_nopiv_percol_trsm_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *tileA, *tileU;
    int                    m, n, k, lda;
    CHAMELEON_Complex64_t *A, *row, pivot;

    tileA = cti_interface_get( descr[0] );
    tileU = cti_interface_get( descr[1] );

    starpu_codelet_unpack_args( cl_arg, &m, &n, &k );

    A   = tileA->mat;
    lda = tileA->ld;
    row = tileU->mat;

    /* Shift to the right column */
    A += k * lda;

    pivot = 1. / *row;
    cblas_zscal( m, CBLAS_SADDR( pivot ), A, 1 );

    /* Update trailing matrix from k+1 to n */
    CORE_zgemm( ChamNoTrans, ChamNoTrans,
                m, n-k-1, 1,
                mzone, A,       lda,
                       row + 1, 1,
                zone,  A + lda, lda );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_nopiv_percol_trsm, cl_zgetrf_nopiv_percol_trsm_cpu_func )

void INSERT_TASK_zgetrf_nopiv_percol_trsm( const RUNTIME_option_t *options,
                                                 int m, int n, int k,
                                                 const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *U, int Um, int Un )
{
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_R(U, Um, Un);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgetrf_nopiv_percol_trsm,
        STARPU_VALUE, &m, sizeof(int),
        STARPU_VALUE, &n, sizeof(int),
        STARPU_VALUE, &k, sizeof(int),

        /* Task handles */
        STARPU_RW, RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,  RTBLKADDR(U, ChamComplexDouble, Um, Un),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgetrf_nopiv_percol_trsm ),
        0 );
}
