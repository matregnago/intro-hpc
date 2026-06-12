/**
 *
 * @file step6.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step6 example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2025-01-29
 *
 */
#include "step6.h"

/*
 * @brief step6 introduces how to use CHAMELEON with MPI.
 * @details This program is a copy of step5 with some additional parameters to
 * be set for the data distribution. To use this program properly CHAMELEON must
 * use StarPU Runtime system and MPI option must be activated at configure.
 * The data distribution used here is 2D block cyclic, see for example
 * http://www.netlib.org/scalapack/slug/node75.html for explanation.
 * The user can enter the parameters of the distribution grid at execution with
 * --p= and --q=
 */
int main(int argc, char *argv[]) {
    size_t N; // matrix order
    int NB;   // number of rows and columns in tiles
    int NRHS; // number of RHS vectors
    int NCPU; // number of cores to use
    int NGPU; // number of gpus (cuda devices) to use
    int GRID_P; // parameter of the 2D block cyclic distribution
    int GRID_Q; // parameter of the 2D block cyclic distribution
    int NMPIPROC = 1; // number of MPI processus
    int UPLO = ChamUpper; // where is stored L

    /* descriptors necessary for calling CHAMELEON tile interface  */
    CHAM_desc_t *descA = NULL, *descAC = NULL, *descB = NULL, *descX = NULL;

    /* declarations to time the program and evaluate performances */
    double fmuls, fadds, flops, gflops, cpu_time;

    /* variable to check the numerical results */
    double anorm, bnorm, xnorm, eps, res;
    int hres;

    /* CHAMELEON sequence uniquely identifies a set of asynchronous function calls
     * sharing common exception handling */
    RUNTIME_sequence_t *sequence = NULL;
    /* CHAMELEON request uniquely identifies each asynchronous function call */
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    /* initialize some parameters with default values */
    int iparam[IPARAM_SIZEOF];
    memset( iparam, 0, sizeof(int) * IPARAM_SIZEOF );
    init_iparam(iparam);

    /* read arguments */
    read_args(argc, argv, iparam);
    N    = iparam[IPARAM_N];
    NB   = iparam[IPARAM_NB];
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
        /* reserve one thread par cuda device to optimize memory transfers */
        iparam[IPARAM_THRDNBR] -= iparam[IPARAM_NCUDAS];
    }
    NCPU = iparam[IPARAM_THRDNBR];
    NGPU = iparam[IPARAM_NCUDAS];

     /* Initialize CHAMELEON with main parameters */
    int rc = CHAMELEON_Init( NCPU, NGPU );
    if (rc != CHAMELEON_SUCCESS) {
        goto finalize;
    }

    /* set some specific parameters related to CHAMELEON: blocks size and inner-blocking size */
    CHAMELEON_Set(CHAMELEON_TILE_SIZE,        iparam[IPARAM_NB] );
    CHAMELEON_Set(CHAMELEON_INNER_BLOCK_SIZE, iparam[IPARAM_IB] );

#if defined(CHAMELEON_USE_MPI)
    NMPIPROC = CHAMELEON_Comm_size();
    /* Check P */
    if ( (iparam[IPARAM_P] > 1) &&
         (NMPIPROC % iparam[IPARAM_P] != 0) ) {
      fprintf(stderr, "ERROR: %d doesn't divide the number of MPI processus %d\n",
              iparam[IPARAM_P], NMPIPROC );
      return EXIT_FAILURE;
    }
#endif
    iparam[IPARAM_Q] = NMPIPROC / iparam[IPARAM_P];
    iparam[IPARAM_NMPI] = NMPIPROC;
    GRID_P = iparam[IPARAM_P];
    GRID_Q = iparam[IPARAM_Q];

    if ( CHAMELEON_Comm_rank() == 0 ){
        /* print informations to user */
        print_header( argv[0], iparam);
    }

    /* Initialize the structure required for CHAMELEON tile interface */
    CHAMELEON_Desc_Create(&descA, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, N, 0, 0, N, N,
                      GRID_P, GRID_Q);
    CHAMELEON_Desc_Create(&descB, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS,
                      GRID_P, GRID_Q);
    CHAMELEON_Desc_Create(&descX, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS,
                      GRID_P, GRID_Q);
    CHAMELEON_Desc_Create(&descAC, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, N, 0, 0, N, N,
                      GRID_P, GRID_Q);

    /* generate A matrix with random values such that it is spd */
    CHAMELEON_dplgsy_Tile( (double)N, ChamUpperLower, descA, 51 );

    /* generate RHS */
    CHAMELEON_dplrnt_Tile( descB, 5673 );

    /* copy A before facto. in order to check the result */
    CHAMELEON_dlacpy_Tile(ChamUpperLower, descA, descAC);

    /* copy B in X before solving
     * same sense as memcpy( X, B, sizeof(double) * N * NRHS ) but for descriptors */
    CHAMELEON_dlacpy_Tile(ChamUpperLower, descB, descX);

    /************************************************************/
    /* solve the system AX = B using the Cholesky factorization */
    /************************************************************/

    cpu_time = -CHAMELEON_timer();

    CHAMELEON_Sequence_Create(&sequence);

    /* Cholesky factorization:
     * A is replaced by its factorization L or L^T depending on uplo */
    CHAMELEON_dpotrf_Tile_Async( UPLO, descA, sequence, &request );

    /* Solve:
     * B is stored in X on entry, X contains the result on exit.
     * Forward and back substitutions
     */
    CHAMELEON_dpotrs_Tile_Async( UPLO, descA, descX, sequence, &request);

    /* Ensure that all data processed on the gpus we are depending on are back
     * in main memory */
    CHAMELEON_Desc_Flush( descA, sequence );
    CHAMELEON_Desc_Flush( descX, sequence );

    /* Synchronization barrier (the runtime ensures that all submitted tasks
     * have been terminated */
    CHAMELEON_Sequence_Wait(sequence);

    rc = sequence->status;
    if ( rc != CHAMELEON_SUCCESS ) {
        fprintf(stderr, "Error in computation (%d)\n", rc);
        goto finalize;
    }
    CHAMELEON_Sequence_Destroy(sequence);

    cpu_time += CHAMELEON_timer();

    /* print informations to user */
    gflops = flops / cpu_time;
    if ( CHAMELEON_Comm_rank() == 0 ) {
        printf( "%9.3f %9.2f\n", cpu_time, gflops);
    }
    fflush( stdout );

    /************************************************************/
    /* check if solve is correct i.e. AX-B = 0                  */
    /************************************************************/

    /* compute norms to check the result */
    anorm = CHAMELEON_dlange_Tile( ChamInfNorm, descAC);
    bnorm = CHAMELEON_dlange_Tile( ChamInfNorm, descB);
    xnorm = CHAMELEON_dlange_Tile( ChamInfNorm, descX);

    /* compute A*X-B, store the result in B */
    CHAMELEON_dgemm_Tile( ChamNoTrans, ChamNoTrans,
                      1.0, descAC, descX, -1.0, descB );
    res = CHAMELEON_dlange_Tile( ChamInfNorm, descB );

    /* check residual and print a message */
    eps = LAPACKE_dlamch_work( 'e' );

    /*
     * if hres = 0 then the test succeed
     * else the test failed
     */
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 );
    if ( CHAMELEON_Comm_rank() == 0 ){
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
    }

    /* deallocate A, B, X, Acpy and associated descriptors descA, ... */
    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descAC );

finalize:
    /*
     * Required semicolon to have at least one inst
     * before the end of OpenMP block.
     */
    ;
    /* Finalize CHAMELEON */
    CHAMELEON_Finalize();

    return EXIT_SUCCESS;
}
