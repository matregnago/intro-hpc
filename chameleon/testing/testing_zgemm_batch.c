/**
 *
 * @file testing_zgemm_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm_batch testing
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
flops_zgemm_batch( int nb, int M, int N, int K )
{
    return flops_zgemm( M, N, K ) * nb;
}

int
testing_zgemm_batch_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async  = parameters_getvalue_int( "async" );
    int                   nb     = run_arg_get_int( args, "nb", 10 );
    int                   ib     = run_arg_get_int( args, "ib", 10 );
    cham_trans_t          transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t          transB = run_arg_get_trans( args, "transB", ChamNoTrans );
    int                   N      = run_arg_get_int( args, "N", 320 );
    int                   M      = run_arg_get_int( args, "M", N );
    int                   K      = run_arg_get_int( args, "K", N );
    CHAMELEON_Complex64_t alpha  = testing_zalea();
    CHAMELEON_Complex64_t beta   = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB  = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC  = run_arg_get_int( args, "seedC", testing_ialea() );

    /* Descriptors */
    int          Am, An, Bm, Bn;
    CHAM_desc_t *descA, *descB, *descC;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta",  beta  );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Calculate the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Am = M;
        An = K;
    }
    else {
        Am = K;
        An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K;
        Bn = N;
    }
    else {
        Bm = N;
        Bn = K;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, Am, An, nb * Am, ib * An, nb * Am, ib * An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, Bm, Bn, nb * Bm, ib * Bn, nb * Bm, ib * Bn );
    parameters_desc_create( "C", &descC, ChamComplexDouble, M,  N,  nb * M,  ib * N,  nb * M,  ib * N  );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Start measurement */
    testing_start( &test_data );
    if ( async ) {
        fprintf( stderr, "Async unavailable yet\n" );
    }
    hres = CHAMELEON_zgemm_batch_Tile( transA, transB, alpha, descA, descB, beta, descC );

    /* Stop measurement */
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm_batch( nb*ib, M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        fprintf( stderr, "Check is not available for zgemm_batch\n" );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );
    parameters_desc_destroy( &descC );

    return hres;
}

testing_t   test_zgemm_batch;
const char *zgemm_batch_params[] = { "nb", "ib", "transA", "transB", "m", "n", "k",
                                     "alpha", "beta", "seedA", "seedB", "seedC", NULL };
const char *zgemm_batch_output[] = { NULL };
const char *zgemm_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgemm_batch_init( void ) __attribute__( ( constructor ) );
void
testing_zgemm_batch_init( void )
{
    test_zgemm_batch.name      = "zgemm_batch";
    test_zgemm_batch.helper    = "Perform nb*ib general matrix-matrix multiply of size MxNxK";
    test_zgemm_batch.params    = zgemm_batch_params;
    test_zgemm_batch.output    = zgemm_batch_output;
    test_zgemm_batch.outchk    = zgemm_batch_outchk;
    test_zgemm_batch.fptr_desc = testing_zgemm_batch_desc;
    test_zgemm_batch.fptr_std  = NULL;
    test_zgemm_batch.next      = NULL;

    testing_register( &test_zgemm_batch );
}
