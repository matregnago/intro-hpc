/**
 *
 * @file step3.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step3 example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @author Lionel Eyraud-Dubois
 * @date 2025-01-29
 *
 */
#include "step3.h"

/*
 * @brief step3 indicates how to give your own tile matrix to CHAMELEON.
 * @details This program is a copy of step2 but instead of using a predefined
 * way for accessing tile data (i.e with CHAMELEON_Desc_Create), we will indicate
 * how to create a CHAMELEON descriptor with an arbitrary tile matrix structure
 * by calling CHAMELEON_Desc_Create_User function.
 * During this step we do not use classical LAPACK matrices (1D array) anymore.
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

    /* Array of pointers to double arrays (representing tiles) */
    double **matA = NULL;

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

    /* allocate tile data */
    matA = allocate_tile_matrix(N, N, NB);

    /*
     * This function is very similar to CHAMELEON_Desc_Create but the way to
     * access matrix tiles can be controled by the user.
     * To do so, three functions with a precise prototype must be given:
     *     - void* get_blkaddr(const CHAM_desc_t *, int, int)
     *     returns a pointer to the tile m, n
     *     - int   get_blkldd (const CHAM_desc_t *, int)
     *     returns the leading dimension of the tile m, n
     *     - int   get_rankof (const CHAM_desc_t *, int, int)
     *     returns the MPI rank of the tile m, n (0 here because we do not
     *     intend to use this program with MPI)
     */
    CHAMELEON_Desc_Create_User( &descA, matA, ChamRealDouble,
                                NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1,
                                user_getaddr_arrayofpointers,
                                user_getblkldd_arrayofpointers,
                                user_getrankof_zero, NULL );

    /*
     * We use the classical CHAMELEON way for accessing tiles for descripotrs
     * B, X and AC. This to show you can define different way to consider tiles
     * in your matrix. The only thing important is to have well defined
     * functions get_blkaddr, get_blkldd, get_rankof corresponding to your data.
     * Note that this call of CHAMELEON_Desc_Create_User routine with
     * chameleon_getaddr_ccrb, chameleon_getblkldd_ccrband chameleon_getrankof_2d functions
     * is equivalent to a call to CHAMELEON_Desc_Create (chameleon_get... are the
     * functions used inside CHAMELEON_Desc_Create).
     */
    CHAMELEON_Desc_Create(&descB, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    CHAMELEON_Desc_Create(&descX, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, 1, 1);
    CHAMELEON_Desc_Create(&descAC, NULL, ChamRealDouble,
                      NB, NB, NB*NB, N, N, 0, 0, N, N, 1, 1);

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

    /* free the matrix of tiles */
    deallocate_tile_matrix(matA, N, N, NB);
    descA->mat = NULL;

    /* free descriptors descA, descB, descX, descAC */
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
