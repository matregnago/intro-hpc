/**
 *
 * @file testing_zpotrf_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf_batch testing
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#endif

static cham_fixdbl_t
flops_zpotrf_batch( int nb, int N )
{
    return flops_zpotrf( N ) * nb;
}

int
testing_zpotrf_batch_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int           async = parameters_getvalue_int( "async" );
    int           nb    = run_arg_get_int( args, "nb", 10 );
    int           ib    = run_arg_get_int( args, "ib", 10 );
    cham_uplo_t   uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int           N     = run_arg_get_int( args, "N", 320 );
    int           seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, N, N, nb * N, ib * N, nb * N, ib * N );

    /* Fill the matrices with random values */
    CHAMELEON_zplghe_batch_Tile( (double)N, descA, seedA );

    /* Start measurement */
    testing_start( &test_data );
    if ( async ) {
        fprintf( stderr, "Async unavailable yet\n" );
    }
    hres = CHAMELEON_zpotrf_batch_Tile( uplo, descA );

    /* Stop measurement */
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf_batch( nb*ib, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        fprintf( stderr, "Check is not available for gemm_batch\n" );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

testing_t   test_zpotrf_batch;
const char *zpotrf_batch_params[] = { "nb", "ib",  "uplo", "n", "seedA", NULL };
const char *zpotrf_batch_output[] = { NULL };
const char *zpotrf_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrf_batch_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrf_batch_init( void )
{
    test_zpotrf_batch.name      = "zpotrf_batch";
    test_zpotrf_batch.helper    = "Perform nb*ib Cholesky factorization potrf( uplo, N, ... )";
    test_zpotrf_batch.params    = zpotrf_batch_params;
    test_zpotrf_batch.output    = zpotrf_batch_output;
    test_zpotrf_batch.outchk    = zpotrf_batch_outchk;
    test_zpotrf_batch.fptr_desc = testing_zpotrf_batch_desc;
    test_zpotrf_batch.fptr_std  = NULL;
    test_zpotrf_batch.next      = NULL;

    testing_register( &test_zpotrf_batch );
}
