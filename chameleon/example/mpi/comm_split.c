/**
 *
 * @file comm_split.c
 *
 * @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon comm_split example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-01-29
 *
 */
#include "comm_split.h"

/*
 * @brief comm_split introduces how to use CHAMELEON with MPI
 * subcommunicators.
 * @details This example shows that Chameleon can be used with custom MPI
 * communicators (different from MPI_COMM_WORLD). Here two different algorithms
 * (potrf and getrf_nopiv) are called at the same time on two different
 * communicators A (0, 2) and B (1, 3). To use this program properly CHAMELEON
 * must use StarPU Runtime system and MPI option must be activated at
 * configure. This program is meant to be run with 4 MPI processes and the data
 * distribution on matrices is 2D block cyclic P=2, Q=1.
 */
int main(int argc, char *argv[]) {

  /* Check that MPI has threads support */
  int thread_support;
  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_support) != MPI_SUCCESS)
  {
    fprintf(stderr,"MPI_Init_thread failed\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if (thread_support == MPI_THREAD_FUNNELED)
  	fprintf(stderr,"Warning: MPI only has funneled thread support, not serialized, hoping this will work\n");
  if (thread_support < MPI_THREAD_FUNNELED)
  	fprintf(stderr,"Warning: MPI does not have thread support!\n");

  /* Check that 4 MPI processes are used */
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if(comm_size != 4)
  {
      printf("This application is meant to be run with 4 MPI processes, not %d.\n", comm_size);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  /* Get my rank in the global communicator */
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Determine the colour and key based on whether my rank is even. */
  char subcommunicator;
  int colour;
  int key = my_rank/2;
  if(my_rank % 2 == 0)
  {
      subcommunicator = 'A';
      colour = 0;
  }
  else
  {
      subcommunicator = 'B';
      colour = 1;
  }

  /* Split the global communicator */
  MPI_Comm new_comm;
  MPI_Comm_split(MPI_COMM_WORLD, colour, key, &new_comm);

  /* Get my rank in the new communicator */
  int my_new_comm_rank;
  MPI_Comm_rank(new_comm, &my_new_comm_rank);

  /* Print my new rank and new communicator */
  printf("[MPI process %d] MPI process %d in subcommunicator %c.\n", my_rank, my_new_comm_rank, subcommunicator);

  /* initialize some parameters with default values */
  int iparam[IPARAM_SIZEOF];
  memset( iparam, 0, sizeof(int) * IPARAM_SIZEOF );
  init_iparam(iparam);

  /* read arguments */
  read_args(argc, argv, iparam);
  int N    = iparam[IPARAM_N]; // matrix order
  int NB   = iparam[IPARAM_NB]; // number of rows and columns in tiles
  int NRHS = iparam[IPARAM_NRHS]; // number of RHS vectors
  int NCPU = iparam[IPARAM_NCPU]; // number of CPU cores to use
  int P = 2;
  int Q = 1;

  /* Initialize CHAMELEON with custom communicator */
  CHAMELEON_InitParComm( NCPU, 0, 1, new_comm );

  CHAMELEON_Set(CHAMELEON_TILE_SIZE, iparam[IPARAM_NB] );

  /* declarations to time the program and evaluate performances */
  double flops, gflops, cpu_time;

  if(my_rank % 2 == 0)
  {
    /* Cholesky on subcommunicator A i.e. 0,2 */

    /* Initialize matrices */
    CHAM_desc_t *descA, *descAC, *descB, *descX;

    CHAMELEON_Desc_Create( &descA, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create( &descAC, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create( &descB, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, P, Q );
    CHAMELEON_Desc_Create( &descX, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, P, Q );

    CHAMELEON_dplgsy_Tile( (double)N, ChamUpperLower, descA, 20 );
    CHAMELEON_dplrnt_Tile( descB, 21 );

    CHAMELEON_dlacpy_Tile( ChamUpperLower, descA, descAC);
    CHAMELEON_dlacpy_Tile( ChamUpperLower, descB, descX);

    cpu_time = -CHAMELEON_timer();

    /* Solve AX = B by Cholesky factorization */
    CHAMELEON_dposv_Tile( ChamLower, descA, descX );

    cpu_time += CHAMELEON_timer();

    flops = flops_dpotrf( N ) + flops_dpotrs( N, NRHS );
    gflops = flops * 1.e-9 / cpu_time;
    if ( CHAMELEON_Comm_rank() == 0 ) {
        printf( "\nCholesky performances on subcommunicator %c:\n", subcommunicator);
        print_header( argv[0], iparam);
        printf( "%9.3f %9.2f\n", cpu_time, gflops);
    }
    fflush( stdout );

    /* compute norms to check the result */
    check( descAC, descX, descB );

    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descAC );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descA );
  }
  else
  {
    /* LU nopiv on subcommunicator B i.e. 1,3 */

    /* Initialize matrices */
    CHAM_desc_t *descA, *descAC, *descB, *descX;

    CHAMELEON_Desc_Create( &descA, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create( &descAC, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create( &descB, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, P, Q );
    CHAMELEON_Desc_Create( &descX, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                           NB, NB, NB*NB, N, NRHS, 0, 0, N, NRHS, P, Q );

    CHAMELEON_dplgtr_Tile( 0        , ChamUpper, descA, 10 );
    CHAMELEON_dplgtr_Tile( (double)N, ChamLower, descA, 11 );
    CHAMELEON_dplrnt_Tile( descB, 12 );

    CHAMELEON_dlacpy_Tile( ChamUpperLower, descA, descAC);
    CHAMELEON_dlacpy_Tile( ChamUpperLower, descB, descX);

    cpu_time = -CHAMELEON_timer();

    /* Solve AX = B by LU factorization */
    CHAMELEON_dgesv_nopiv_Tile( descAC, descX );

    cpu_time += CHAMELEON_timer();

    flops = flops_dgetrf( N, N ) + flops_dgetrs( N, NRHS );
    gflops = flops * 1.e-9 / cpu_time;
    if ( CHAMELEON_Comm_rank() == 0 ) {
        printf( "\nLU nopiv performances on subcommunicator %c:\n", subcommunicator);
        print_header( argv[0], iparam);
        printf( "%9.3f %9.2f\n", cpu_time, gflops);
    }
    fflush( stdout );

    /* compute norms to check the result */
    check( descA, descX, descB );

    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descAC );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descA );
  }

  CHAMELEON_Finalize();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
