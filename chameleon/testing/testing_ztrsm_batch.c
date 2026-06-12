/**
 *
 * @file testing_ztrsm_batch.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm_batch testing
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
flops_ztrsm_batch( int nb, cham_side_t side, int M, int N )
{
    return flops_ztrsm( side, M, N ) * nb;
}

int
testing_ztrsm_batch_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_int( args, "nb", 10 );
    int                   ib    = run_arg_get_int( args, "ib", 10 );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N     = run_arg_get_int( args, "N", 320 );
    int                   M     = run_arg_get_int( args, "M", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB;

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Calculate the dimensions according to the transposition */
    if ( side == ChamLeft ) {
        Am = M;
        An = M;
    }
    else {
        Am = N;
        An = N;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, Am, An, nb * Am, ib * An, nb * Am, ib * An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, M,  N,  nb * M,  ib * N,  nb * M,  ib * N  );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Start measurement */
    testing_start( &test_data );
    if ( async ) {
        fprintf( stderr, "Async unavailable yet\n" );
    }
    hres = CHAMELEON_ztrsm_batch_Tile( side, uplo, trans, diag, alpha, descA, descB );

     /* Stop measurement */
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrsm_batch( nb*ib, side, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        fprintf( stderr, "Check is not available for ztrsm_batch\n" );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}

testing_t   test_ztrsm_batch;
const char *ztrsm_batch_params[] = { "nb", "ib", "side",  "uplo", "trans", "diag", "m", "n",
                                     "alpha", "seedA", "seedB", NULL };
const char *ztrsm_batch_output[] = { NULL };
const char *ztrsm_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrsm_batch_init( void ) __attribute__( ( constructor ) );
void
testing_ztrsm_batch_init( void )
{
    test_ztrsm_batch.name      = "ztrsm_batch";
    test_ztrsm_batch.helper    = "Perform nb*ib triangular solve trsm( side, uplo, trns, diag, M, N, ... )";
    test_ztrsm_batch.params    = ztrsm_batch_params;
    test_ztrsm_batch.output    = ztrsm_batch_output;
    test_ztrsm_batch.outchk    = ztrsm_batch_outchk;
    test_ztrsm_batch.fptr_desc = testing_ztrsm_batch_desc;
    test_ztrsm_batch.fptr_std  = NULL;
    test_ztrsm_batch.next      = NULL;

    testing_register( &test_ztrsm_batch );
}
