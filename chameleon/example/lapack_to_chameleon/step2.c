/**
 *
 * @file step2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step2 example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2025-01-29
 *
 */
#include "step2.h"

/*
 * @brief step2 introduces the CHAMELEON tile interface.
 * @details This program is a copy of step1 but instead of using the LAPACK
 * interface which leads to copy LAPACK matrices inside CHAMELEON routines we use
 * the tile interface.
 * We will still use standard format of matrix but we will see how to give this
 * matrix to create a CHAMELEON descriptor, a structure wrapping data on which we
 * want to apply sequential task based algorithms.
 * The goal is also to show examples of calling CHAMELEON tile interface.
 */
int main(int argc, char *argv[]) {
    size_t N; // matrix order
    int NB;   // number of rows and columns in tiles
    int NRHS; // number of RHS vectors
    int NCPU; // number of cores to use
    int NGPU; // number of gpus (cuda devices) to use
    int UPLO = ChamUpper; // where is stored L

    /* descriptors necessary for calling CHAMELEON tile interface  */
    CHAM_desc_t *descA = NULL, *descAC = NULL, *descB = NULL, *descX = NULL;

    /* declarations to time the program and evaluate performances */
    double fmuls, fadds, flops, gflops, cpu_time;

    /* variable to check the numerical results */
    double anorm, bnorm, xnorm, eps, res;
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
    NCPU = iparam[IPARAM_THRDNBR];
    NGPU = 0;

    /* print informations to user */
    print_header( argv[0], iparam);

    /* Initialize CHAMELEON with main parameters */
    int rc = CHAMELEON_Init( NCPU, NGPU );
    if (rc != CHAMELEON_SUCCESS) {
        goto finalize;
    }

    /* Question chameleon to get the block (tile) size (number of columns) */
    CHAMELEON_Get( CHAMELEON_TILE_SIZE, &NB );

    /*
     * Allocate memory for our data using a C macro (see step2.h)
     *     - matrix A                   : size N x N
     *     - set of RHS vectors B       : size N x NRHS
     *     - set of solutions vectors X : size N x NRHS
     */
    double *A    = malloc( sizeof(double) * N * N    );
    double *Acpy = malloc( sizeof(double) * N * N    );
    double *B    = malloc( sizeof(double) * N * NRHS );
    double *X    = malloc( sizeof(double) * N * NRHS );

    /*
     * Initialize the structure required for CHAMELEON tile interface
     * CHAM_desc_t is a structure wrapping your data allowing CHAMELEON to get
     * pointers to tiles. A tile is a data subset of your matrix on which we
     * apply some optimized CPU/GPU kernels.
     * Notice that this routine suppose your matrix is a contiguous vector of
     * data (1D array), as a data you would give to BLAS/LAPACK.
     * Main arguments:
     *     - descA is a pointer to a descriptor, you need to give the address
     *     of this pointer
     *     - if you want to give your allocated matrix give its address,
     *     if not give a NULL pointer, the routine will allocate the memory
     *     and you access the matrix data with descA->mat
     *     - give the data type (ChamByte, ChamInteger, ChamRealFloat,
     *     ChamRealDouble, ChamComplexFloat, ChamComplexDouble)
     *     - number of rows in a block (tile)
     *     - number of columns in a block (tile)
     *     - number of elements in a block (tile)
     * The other parameters are specific, use:
     * CHAMELEON_Desc_Create( ... , 0, 0, number of rows, number of columns, 1, 1);
     * Have a look to the documentation for details about these parameters.
     */
    CHAMELEON_Desc_Create(&descA, A, ChamRealDouble,
                      NB, NB,  NB*NB, N, N, 0, 0, N, N, 1, 1);
    CHAMELEON_Desc_Create(&descB, B, ChamRealDouble,
                      NB, NB,  NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    CHAMELEON_Desc_Create(&descX, X, ChamRealDouble,
                      NB, NB,  NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    CHAMELEON_Desc_Create(&descAC, Acpy, ChamRealDouble,
                      NB, NB,  NB*NB, N, N, 0, 0, N, N, 1, 1);

    /* You could alternatively create descriptors wrapping your allocated
     * matrices to avoid copies Lap2Desc with the following */
    //CHAMELEON_Desc_Create(&descA, A, ChamRealDouble,
    //                  NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1);
    /* Be aware that for distributed data (with MPI) you need to respect a
     * precise partitionning called 2D cyclic (see the documentation) */

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

    /* Cholesky factorization:
     * A is replaced by its factorization L or L^T depending on uplo */
    CHAMELEON_dpotrf_Tile( UPLO, descA );

    /* Solve:
     * B is stored in X on entry, X contains the result on exit.
     * Forward and back substitutions
     */
    CHAMELEON_dpotrs_Tile( UPLO, descA, descX );

    cpu_time += CHAMELEON_timer();

    /* print informations to user */
    gflops = flops / cpu_time;
    printf( "%9.3f %9.2f\n", cpu_time, gflops);
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

    return rc;
}
