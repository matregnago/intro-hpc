/**
 *
 * @file testing_zgram.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram testing
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgram( int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t N     = _N;

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* 5 multiplications and 3 addition per element */
    flops = ( 5. * 6. + 6. ) * N * N;
#else
    flops = ( 5. + 3. ) * N * N;
#endif

    return flops;
}

int
testing_zgram_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Fill the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    if ( async ) {
        ws = CHAMELEON_zgram_WS_Alloc( descA );
    }

    /* Compute the gram matrix transformation */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgram_Tile_Async( uplo, descA, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgram_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgram( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    parameters_desc_destroy( &descA );

    (void)check;
    return hres;
}

int
testing_zgram_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fill the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N,A, LDA, seedA );

    /* Compute the gram matrix transformation */
    testing_start( &test_data );
    hres = CHAMELEON_zgram( uplo, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgram( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zgram;
const char *zgram_params[] = { "mtxfmt", "nb", "uplo", "n", "n", "lda", "seedA", NULL };
const char *zgram_output[] = { NULL };
const char *zgram_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgram_init( void ) __attribute__( ( constructor ) );
void
testing_zgram_init( void )
{
    test_zgram.name   = "zgram";
    test_zgram.helper = "General Gram matrix transformation";
    test_zgram.params = zgram_params;
    test_zgram.output = zgram_output;
    test_zgram.outchk = zgram_outchk;
    test_zgram.fptr_desc = testing_zgram_desc;
    test_zgram.fptr_std  = testing_zgram_std;
    test_zgram.next   = NULL;

    testing_register( &test_zgram );
}
