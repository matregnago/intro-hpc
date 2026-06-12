/**
 *
 * @file comm_split.h
 *
 * @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon comm_split example header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-05-06
 *
 */
#ifndef _comm_split_h_
#define _comm_split_h_

#include <coreblas/lapacke.h>
#include <chameleon.h>
#include <chameleon/timer.h>
#include <chameleon/flops.h>
#include <mpi.h>
#include <string.h>

/* Integer parameters for comm_split */
enum iparam_comm_split {
    IPARAM_NCPU, /* Number of CPUs                       */
    IPARAM_N,    /* Number of columns/rows of the matrix */
    IPARAM_NB,   /* Number of columns/rows in a tile     */
    IPARAM_NRHS, /* Number of RHS                        */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in comm_split.c main program */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_NCPU ] = -1;
    iparam[IPARAM_N    ] = 1000;
    iparam[IPARAM_NB   ] = 500;
    iparam[IPARAM_NRHS ] = 1;
 }

/**
 * Print how to use the program
 */
static void show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help   Show this help\n"
            "\n"
            "  --n=X    dimension (N). (default: 1000)\n"
            "  --nb=X   NB size. (default: 500)\n"
            "  --nrhs=X number of RHS. (default: 1)\n"
            "\n"
            "  --cpus=X Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "\n");
}

static int startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

/**
 * Read arguments following comm_split program call
 */
static void read_args(int argc, char *argv[], int *iparam){
    int i;
    for (i = 1; i < argc && argv[i]; ++i) {
        if ( startswith( argv[i], "--help") || startswith( argv[i], "-help") ||
             startswith( argv[i], "--h") || startswith( argv[i], "-h") ) {
            show_help( argv[0] );
            exit(0);
        } else if (startswith( argv[i], "--n=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_N]) );
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NB]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--cpus=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NCPU]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }
}

/**
 * Print a header message to summarize main parameters
 */
static void print_header(char *prog_name, int * iparam) {
    double    eps = LAPACKE_dlamch_work( 'e' );
    printf( "#\n"
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb cpu: %d\n"
            "# N:      %d\n"
            "# NB:     %d\n"
            "# eps:    %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_NCPU],
            iparam[IPARAM_N],
            iparam[IPARAM_NB],
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

/**
 * Check that AX=B is correct
 */
static void check(CHAM_desc_t *A, CHAM_desc_t *X, CHAM_desc_t *B) {
  double anorm = CHAMELEON_dlange_Tile( ChamInfNorm, A);
  double xnorm = CHAMELEON_dlange_Tile( ChamInfNorm, X);
  double bnorm = CHAMELEON_dlange_Tile( ChamInfNorm, B);
  CHAMELEON_dgemm_Tile( ChamNoTrans, ChamNoTrans,
                        1.0, A, X, -1.0, B );
  double res = CHAMELEON_dlange_Tile( ChamInfNorm, B );
  double eps = LAPACKE_dlamch_work( 'e' );
  int N = X->lm;
  int hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 );
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
  return;
}

#endif /* _comm_split_h_ */
