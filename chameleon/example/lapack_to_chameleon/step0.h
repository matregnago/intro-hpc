/**
 *
 * @file step0.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step0 example header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-09-19
 *
 */
#ifndef _step0_h_
#define _step0_h_

/* Common include for all steps of the tutorial */
#include "lapack_to_chameleon.h"

/* Specific includes for step 0 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>

/* Integer parameters for step0 */
enum iparam_step0 {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_LDA,            /* Leading dimension of A                     */
    IPARAM_LDB,            /* Leading dimension of B                     */
    IPARAM_IB,             /* Inner-blocking size                        */
    IPARAM_NRHS,           /* Number of RHS                              */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in step0.c main program */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_LDA           ] = iparam[IPARAM_N];
    iparam[IPARAM_LDB           ] = iparam[IPARAM_N];
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NRHS          ] = 1;
 }

/**
 * Print how to use the program
 */
static void show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --n=X            dimension (N) .(default: 500)\n"
            "  --ib=X           IB size. (default: 32)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "\n");
}

/**
 * Read arguments following step0 program call
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
        } else if (startswith( argv[i], "--ib=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_IB]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }
}

/**
 * Print a header message to summarize main parameters
 */
static void print_header(char *prog_name, int * iparam) {
    double eps = LAPACKE_dlamch_work( 'e' );

    printf( "#\n"
            "# CHAMELEON %s\n"
            "# Nb threads: %d\n"
            "# N:          %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            prog_name,
            iparam[IPARAM_THRDNBR],
            iparam[IPARAM_N],
            iparam[IPARAM_IB],
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

/* The following is a copy paste from coreblas/core_zplgsy and core_zplrnt.
 * We need these routines to generate the matrices but we don't want to link with chameleon here.
 */
#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20
#define NBELEM   1

static unsigned long long int
Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, i++) {
    if (n & 1)
      ran = a_k * ran + c_k;
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

/**
 * CORE_dplgsy - Generate a tile for random symmetric (positive definite if 'bump' is large enough) matrix.
 */
static void CORE_dplgsy( double bump, int m, int n, double *A, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    double *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)bigM;

    /*
     * Tile diagonal
     */
    if ( m0 == n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = j; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
                tmp++;
            }
            tmp  += (lda - i + j + 1);
            jump += bigM + 1;
        }

        for (j = 0; j < n; j++) {
            A[j+j*lda] += bump;

            for (i=0; i<j; i++) {
                A[lda*j+i] = A[lda*i+j];
            }
        }
    }
    /*
     * Lower part
     */
    else if ( m0 > n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = 0; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
                tmp++;
            }
            tmp  += (lda - i);
            jump += bigM;
        }
    }
    /*
     * Upper part
     */
    else if ( m0 < n0 ) {
        /* Overwrite jump */
        jump = (unsigned long long int)n0 + (unsigned long long int)m0 * (unsigned long long int)bigM;

        for (i = 0; i < m; i++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (j = 0; j < n; j++) {
                A[j*lda+i] = 0.5f - ran * RndF_Mul;
            }
            jump += bigM;
        }
    }
}

/**
 * CORE_dplrnt - Generate a tile for random matrix.
 */
static void CORE_dplrnt( int m, int n, double *A, int lda,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
	double *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)bigM;

    for (j=0; j<n; ++j ) {
        ran = Rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            *tmp = 0.5f - ran * RndF_Mul;
            ran  = Rnd64_A * ran + Rnd64_C;
            tmp++;
        }
        tmp  += lda-i;
        jump += bigM;
    }
}

#endif /* _step0_h_ */
