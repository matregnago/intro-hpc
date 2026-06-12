/**
 *
 * @file step3.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon step3 example header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Pierre Esterie
 * @date 2025-01-29
 *
 */
#ifndef _step3_h_
#define _step3_h_

/* Common include for all steps of the tutorial */
#include "lapack_to_chameleon.h"

/* Specific includes for step 3 */
#include <coreblas/lapacke.h>
#include <chameleon.h>

/* Integer parameters for step3 */
enum iparam_step3 {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_NRHS,           /* Number of RHS                              */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in step3.c main program */

/**
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_N             ] = 500;
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
            "  --n=X            dimension (N). (default: 500)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "\n");
}

/**
 * Read arguments following step3 program call
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
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
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
            0,
            iparam[IPARAM_N],
            128,
            32,
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

/**
 *  Function that allocate an array of pointers to square tiles (allocated to 0)
 */
double **allocate_tile_matrix(int m, int n, int nb){
    int i;
    int mt, nt;
    double **mat;

    /* compute number of tiles in rows and columns */
    mt = (m%nb==0) ? (m/nb) : (m/nb+1);
    nt = (n%nb==0) ? (n/nb) : (n/nb+1);
    mat = malloc( sizeof(double*) * mt * nt );
    if (!mat){
        printf ("\nIn allocate_tile_matrix, memory Allocation Failure of mat !\n\n");
        exit (EXIT_FAILURE);
    }
    for (i=0; i<mt*nt; i++){
        *(mat+i) = calloc( nb*nb, sizeof(double) );
        if (!*(mat+i)){
            printf ("\nIn allocate_tile_matrix, memory Allocation Failure of *(mat+i) !\n\n");
            exit (EXIT_FAILURE);
        }
    }
    return mat;
}

/**
 *  Function that deallocate an array of pointers to square tiles
 */
static void deallocate_tile_matrix(double **mat, int m, int n, int nb){
    int i;
    int mt, nt;

    /* compute number of tiles in rows and columns */
    mt = (m%nb==0) ? (m/nb) : (m/nb+1);
    nt = (n%nb==0) ? (n/nb) : (n/nb+1);
    for (i=0; i<mt*nt; i++) free(*(mat+i));
    free(mat);
}

/**
 *  Function to return address of block (m,n)
 */
inline static void* user_getaddr_arrayofpointers(const CHAM_desc_t *A, int m, int n)
{
    double **matA = (double **)A->mat;
    size_t mm = (size_t)m + (size_t)A->i / A->mb;
    size_t nn = (size_t)n + (size_t)A->j / A->nb;
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / chameleon_desc_datadist_get_iparam(A, 0);
    nn = nn / chameleon_desc_datadist_get_iparam(A, 1);
#endif

    offset = A->mt*nn + mm;
    return (void*)( *(matA + offset) );
}

/**
 *  Function to return the leading dimension of element A(m,*)
 */
inline static int user_getblkldd_arrayofpointers(const CHAM_desc_t *A, int m)
{
    (void)m;
    return A->mb;
}

/**
 *  Function to return MPI rank of element A(m,n)
 */
inline static int user_getrankof_zero(const CHAM_desc_t *A, int m, int n)
{
    (void)A; (void)m; (void)n;
    return 0;
}

#endif /* _step3_h_ */
