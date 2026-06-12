/**
 *
 * @file testing_zplrnk.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnk testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zplrnk_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int async = parameters_getvalue_int( "async" );
    int nb    = run_arg_get_nb( args );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int K     = run_arg_get_int( args, "K", N );
    int LDC   = run_arg_get_int( args, "LDC", M );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descC;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    parameters_desc_create( "C", &descC, ChamComplexDouble, nb, nb, LDC, N, M, N );

    /* Calculates the random rank-k matrix */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zplrnk_Tile_Async( K, descC, seedA, seedB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zplrnk_Tile( K, descC, seedA, seedB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_zrankk( args, K, descC );
    }

    parameters_desc_destroy( &descC );

    return hres;
}

int
testing_zplrnk_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int nb    = run_arg_get_nb( args );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int K     = run_arg_get_int( args, "K", N );
    int LDC   = run_arg_get_int( args, "LDC", M );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *C;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    C = malloc( sizeof(CHAMELEON_Complex64_t) * LDC * N );

    /* Calculates the random rank-k matrix */
    testing_start( &test_data );
    hres = CHAMELEON_zplrnk( M, N, K, C, LDC, seedA, seedB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_zrankk_std( args, M, N, K, C, LDC );
    }

    free( C );

    return hres;
}

testing_t   test_zplrnk;
const char *zplrnk_params[] = { "nb", "m", "n", "k", "ldc", "seedA", "seedB", NULL };
const char *zplrnk_output[] = { NULL };
const char *zplrnk_outchk[] = { "||A||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zplrnk_init( void ) __attribute__( ( constructor ) );
void
testing_zplrnk_init( void )
{
    test_zplrnk.name   = "zplrnk";
    test_zplrnk.helper = "General rank-k matrix generation";
    test_zplrnk.params = zplrnk_params;
    test_zplrnk.output = zplrnk_output;
    test_zplrnk.outchk = zplrnk_outchk;
    test_zplrnk.fptr_desc = testing_zplrnk_desc;
    test_zplrnk.fptr_std  = testing_zplrnk_std;
    test_zplrnk.next   = NULL;

    testing_register( &test_zplrnk );
}
