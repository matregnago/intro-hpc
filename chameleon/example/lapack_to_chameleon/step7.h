/**
 *
 * @file step7.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step7 example header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @date 2025-01-24
 *
 */
#ifndef _step7_h_
#define _step7_h_

/* Common include for all steps of the tutorial */
#include "lapack_to_chameleon.h"

/* Specific includes for step 7 */
#include <coreblas/lapacke.h>
#include <coreblas.h>
#include <chameleon.h>
// #if defined(CHAMELEON_USE_MPI)
// #include <mpi.h>
// #endif

/* Integer parameters for step7 */
enum iparam_step7 {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_NCUDAS,         /* Number of cuda devices                     */
    IPARAM_NMPI,           /* Number of cuda devices                     */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_NB,             /* Number of columns in a tile                */
    IPARAM_IB,             /* Inner-blocking size                        */
    IPARAM_NRHS,           /* Number of RHS                              */
    IPARAM_P,              /* 2D block cyclic distribution parameter MB  */
    IPARAM_Q,              /* 2D block cyclic distribution parameter NB  */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in step7.c main program */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_NCUDAS        ] = 0;
    iparam[IPARAM_NMPI          ] = 1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_NB            ] = 128;
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NRHS          ] = 1;
    iparam[IPARAM_P             ] = 1;
    iparam[IPARAM_Q             ] = 1;
 }

/**
 * Function used to build matrix blocks
 * Cham_build_plgsy_cpu : random symmetric positive definite
 * Cham_build_plrnt_cpu : random general
 * These 2 functions use data_pl to get data on the matrix to build, passed through the opaque pointer 'op_args'
 * The function is expected to build the block of matrix (m,n) and store it in tileA
 */
struct data_pl {
  double                 bump;
  unsigned long long int seed;
};

static int Cham_build_plgsy_cpu( void *op_args, cham_uplo_t uplo, int m, int n, int ndata,
                                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    struct data_pl *data = (struct data_pl *)op_args;
    int             tempmm, tempnn;

    /* Get the dimension of the tile */
    tempmm = descA->get_blkdim( descA, m, DIM_m, descA->m );
    tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );

    /* fill the tile with the coreblas function plgsy = random SPD matrix generator */
    TCORE_dplgsy( data->bump, tempmm, tempnn, tileA,
                  descA->m, m * descA->mb, n * descA->nb, data->seed );

    (void)uplo;
    (void)ndata;
    return 0;
}

static int Cham_build_plrnt_cpu( void *op_args, cham_uplo_t uplo, int m, int n, int ndata,
                                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    struct data_pl *data = (struct data_pl *)op_args;
    int             tempmm, tempnn;

    /* Get the dimension of the tile */
    tempmm = descA->get_blkdim( descA, m, DIM_m, descA->m );
    tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );

    /* fill the tile with the coreblas function plrnt = random general matrix generator */
    TCORE_dplrnt( tempmm, tempnn, tileA,
                  descA->m, m * descA->mb, n * descA->nb, data->seed );

    (void)uplo;
    (void)ndata;
    return 0;
}

/**
 * Print how to use the program
 */
static void show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --n=X            dimension (N). (default: 500)\n"
            "  --nb=X           NB size. (default: 128)\n"
            "  --ib=X           IB size. (default: 32)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "  --p=X            2D block cyclic distribution parameter MB. (default: 1)\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "  --gpus=X         Number of GPU workers (default: 0)\n"
            "\n");
}

/**
 * Read arguments following step7 program call
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
        } else if (startswith( argv[i], "--ib=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_IB]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--p=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_P]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else if (startswith( argv[i], "--gpus=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NCUDAS]) );
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
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb gpus:    %d\n"
            "# Nb mpi:     %d\n"
            "# PxQ:        %dx%d\n"
            "# N:          %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            iparam[IPARAM_NCUDAS],
            iparam[IPARAM_NMPI],
            iparam[IPARAM_P], iparam[IPARAM_Q],
            iparam[IPARAM_N],
            iparam[IPARAM_NB],
            iparam[IPARAM_IB],
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

#endif /* _step7_h_ */
