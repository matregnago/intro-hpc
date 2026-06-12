/**
 *
 * @file step0.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step0 example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include "step0.h"

/*
 * @brief step0 is a simple Cholesky example using the C interface of
 * BLAS/LAPACK.
 * @details This program solves a linear system AX=B with matrix A symmetric
 * positive definite.
 * The matrix A is first factorized using the Cholesky factorization, A = LL^T.
 * Then the solution X is calculated thanks to forward and back substitutions.
 * We use the C interface of BLAS and LAPACK, that is, CBLAS and LAPACKE, to solve this problem.
 * This program is either sequential if it is linked with a sequential LAPACK library
 * or parallel if it is linked with a multithreaded LAPACK library.
 * The precision is: double
 */
int main(int argc, char *argv[]) {
    size_t N;    // matrix order
    size_t NRHS; // number of RHS vectors

    /* declarations to time the program and evaluate performances */
    double fmuls, fadds, flops, gflops, cpu_time, eps, res;
    int hres;

    /* initialize some parameters with default values */
    int iparam[IPARAM_SIZEOF];
    memset( iparam, 0, sizeof(int) * IPARAM_SIZEOF );
    init_iparam(iparam);

    /* read arguments */
    read_args(argc, argv, iparam);
    N    = iparam[IPARAM_N];
    NRHS = iparam[IPARAM_NRHS];

    /* compute the algorithm complexity to evaluate performances */
    fadds = (double)( FADDS_POTRF(N) + 2 * FADDS_TRSM(N,NRHS) );
    fmuls = (double)( FMULS_POTRF(N) + 2 * FMULS_TRSM(N,NRHS) );
    flops = 1e-9 * (fmuls + fadds);

    /* initialize the number of thread if not given by the user in argv
     * It makes sense only if this program is linked with pthread and
     * multithreaded BLAS and LAPACK */
    if ( iparam[IPARAM_THRDNBR] == -1 ) {
        get_thread_count( &(iparam[IPARAM_THRDNBR]) );
    }

    /* print informations to user */
    print_header( argv[0], iparam);

    /*
     * allocate memory for our data using a C macro (see step0.h)
     *     - matrix A                   : size N x N
     *     - set of RHS vectors B       : size N x NRHS
     *     - set of solutions vectors X : size N x NRHS
     */
    double *A    = malloc( sizeof(double) * N * N    );
    double *Acpy = malloc( sizeof(double) * N * N    );
    double *B    = malloc( sizeof(double) * N * NRHS );
    double *X    = malloc( sizeof(double) * N * NRHS );

    /* generate A matrix with random values such that it is spd */
    CORE_dplgsy( (double)N, N, N, A, N, N, N, N, 51 );

    /* generate RHS */
    CORE_dplrnt( N, NRHS, B, N, N, N, NRHS, 5673 );

    /* copy A before facto. in order to check the result */
    memcpy( Acpy, A, sizeof(double) * N * N );

    /* copy B in X before solving */
    memcpy( X, B, sizeof(double) * N * NRHS );

    /************************************************************/
    /* solve the system AX = B using the Cholesky factorization */
    /************************************************************/

    cpu_time = -CHAMELEON_timer();

    /* Cholesky factorization:
     * A is replaced by its factorization L or L^T depending on uplo */
    hres = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', N, A, N );
    if ( hres != 0 ) {
        fprintf( stderr, "Error during Cholesky factorization call (%d)\n", hres );
        goto end;
    }

    /* Solve:
     * B is stored in X on entry, X contains the result on exit.
     * Forward ...
     */
    cblas_dtrsm(
        CblasColMajor,
        CblasLeft,
        CblasUpper,
        CblasConjTrans,
        CblasNonUnit,
        N, NRHS, 1.0, A, N, X, N);
    /* ... and back substitution */
    cblas_dtrsm(
        CblasColMajor,
        CblasLeft,
        CblasUpper,
        CblasNoTrans,
        CblasNonUnit,
        N, NRHS, 1.0, A, N, X, N);

    cpu_time += CHAMELEON_timer();

    /* print informations to user */
    gflops = flops / cpu_time;
    printf( "%9.3f %9.2f\n", cpu_time, gflops);
    fflush( stdout );

    /************************************************************/
    /* check if solve is correct i.e. AX-B = 0                  */
    /************************************************************/

    /* compute norms to check the result */
    double anorm = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', N, N, Acpy, N);
    double bnorm = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', N, NRHS, B, N);
    double xnorm = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', N, NRHS, X, N);

    /* compute A*X-B, store the result in B */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                N, NRHS, N, 1.0, Acpy, N, X, N, -1.0, B, N);
    res = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', N, NRHS, B, N);

    /* check residual and print a message */
    eps = LAPACKE_dlamch_work( 'e' );

    /*
     * if hres = 0 then the test succeed
     * else the test failed
     */
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 );
    printf( "   ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN\n");
    if (hres) {
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e FAILURE \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));
    }
    else {
        printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e SUCCESS \n",
            res, anorm, xnorm, bnorm,
            res / N / eps / (anorm * xnorm + bnorm ));
    }

  end:
    /* deallocate data */
    free(A);
    free(Acpy);
    free(B);
    free(X);

    return EXIT_SUCCESS;
}
