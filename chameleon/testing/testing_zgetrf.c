/**
 *
 * @file testing_zgetrf.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @author Xavier Lacoste
 * @author Florent Pruvost
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <chameleon/getenv.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/lapacke.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zgetrf_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;
    int        P, Q;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    int         ib    = run_arg_get_ib( args );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         minMN = chameleon_min( M, N );

    /* Descriptors */
    CHAM_desc_t *descA;
    CHAM_ipiv_t *descIPIV;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );
    if ( descA->get_rankof_init != chameleon_getrankof_2d ) {
        parameters_desc_destroy( &descA );
        fprintf( stderr,
                 "SKIPPED: zgetrf kernel only supports 2DBC data distributions\n" );
        return -1;
    }

    P = chameleon_desc_datadist_get_iparam( descA, 0 );
    Q = chameleon_desc_datadist_get_iparam( descA, 1 );

    CHAMELEON_Ipiv_Create( &descIPIV, ChamLeft, descA->mb, minMN, P, P*Q, NULL );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    if ( async ) {
        ws = CHAMELEON_zgetrf_WS_Alloc( descA );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgetrf_Tile_Async( descA, descIPIV, ws, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgetrf_Tile( descA, descIPIV );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorization and residual */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zlaswp_Tile( ChamLeft, ChamDirForward, descA0, 1, minMN, descIPIV );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    if ( ws != NULL ) {
        CHAMELEON_zgetrf_WS_Free( ws );
    }

    CHAMELEON_Ipiv_Destroy( &descIPIV );
    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zgetrf_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int         api   = parameters_getvalue_int( "api" );
#endif
    int         nb    = run_arg_get_nb( args );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         minMN = chameleon_min( M, N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;
    int *IPIV;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    IPIV = malloc( sizeof(int) * minMN );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zgetrf( LAPACK_COL_MAJOR, M, N, A, LDA, IPIV );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zgetrf( M, N, A, LDA, IPIV );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zgetrf( CblasColMajor, M, N, A, LDA, IPIV );
        break;
#endif
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

#if !defined(CHAMELEON_SIMULATION)
    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplrnt( M, N, A0, LDA, seedA );

        /* Compute the permutation of A0: P * A0 */
        CHAMELEON_zlaswp( ChamLeft, ChamDirForward, M, N, A0, LDA, 1, minMN, IPIV );

        hres += check_zxxtrf_std( args, ChamGeneral, ChamUpperLower, M, N, A0, A, LDA );

        free( A0 );
    }
#endif
#endif

    free ( IPIV );
    free( A );

    (void)check;
    return hres;
}

testing_t   test_zgetrf;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zgetrf_params[] = { "m", "n", "lda", "seedA", NULL };
#else
const char *zgetrf_params[] = { "mtxfmt", "nb", "ib", "m", "n", "lda", "seedA", NULL };
#endif
const char *zgetrf_output[] = { NULL };
const char *zgetrf_outchk[] = { "||A||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrf_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrf_init( void )
{
    test_zgetrf.name   = "zgetrf";
    test_zgetrf.helper = "General LU factorization (with partial pivoting)";
    test_zgetrf.params = zgetrf_params;
    test_zgetrf.output = zgetrf_output;
    test_zgetrf.outchk = zgetrf_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zgetrf.fptr_desc = NULL;
#else
    test_zgetrf.fptr_desc = testing_zgetrf_desc;
#endif
    test_zgetrf.fptr_std  = testing_zgetrf_std;
    test_zgetrf.next   = NULL;

    testing_register( &test_zgetrf );
}
