/**
 *
 * @file testing_zcesca.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca testing
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
flops_zcesca( int _M, int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t M     = _M;
    cham_fixdbl_t N     = _N;

#if defined(PRECISION_z) || defined(PRECISION_c)
    /*  2 multiplications and 5 addition per element */
    flops = ( 2. * 6. + 10. ) * M * N;
#else
    flops = ( 2. + 5. ) * M * N;
#endif

    return flops;
}

int
testing_zcesca_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int async = parameters_getvalue_int( "async" );
    int nb    = run_arg_get_nb( args );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int LDA   = run_arg_get_int( args, "LDA", M );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    if ( async ) {
        ws = CHAMELEON_zcesca_WS_Alloc( descA );
    }

    /* Compute the centered-scaled matrix transformation */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zcesca_Tile_Async( 1, 1, ChamColumnwise, descA, ws,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zcesca_Tile( 1, 1, ChamColumnwise, descA, NULL, NULL );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zcesca( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    parameters_desc_destroy( &descA );

    (void)check;
    return hres;
}

int
testing_zcesca_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int nb    = run_arg_get_nb( args );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int LDA   = run_arg_get_int( args, "LDA", M );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA * N );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Compute the centered-scaled matrix transformation */
    testing_start( &test_data );
    hres = CHAMELEON_zcesca( 1, 1, ChamColumnwise, M, N, A, LDA, NULL, NULL );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zcesca( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zcesca;
const char *zcesca_params[] = { "mtxfmt", "nb", "trans", "m", "n", "lda", "seedA", NULL };
const char *zcesca_output[] = { NULL };
const char *zcesca_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zcesca_init( void ) __attribute__( ( constructor ) );
void
testing_zcesca_init( void )
{
    test_zcesca.name   = "zcesca";
    test_zcesca.helper = "General centered-scaled matrix transformation";
    test_zcesca.params = zcesca_params;
    test_zcesca.output = zcesca_output;
    test_zcesca.outchk = zcesca_outchk;
    test_zcesca.fptr_desc = testing_zcesca_desc;
    test_zcesca.fptr_std  = testing_zcesca_std;
    test_zcesca.next   = NULL;

    testing_register( &test_zcesca );
}
