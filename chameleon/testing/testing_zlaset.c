/**
 *
 * @file testing_zlaset.c
 *
 * @copyright 2019-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset testing
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zlaset_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpperLower );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();

    /* Descriptors */
    CHAM_desc_t *descA;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );

    /* Apply laset */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlaset_Tile_Async( uplo, alpha, beta, descA,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlaset_Tile( uplo, alpha, beta, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, 0.0 );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zset( args, uplo, alpha, beta, descA );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zlaset_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   nb    = run_arg_get_nb( args );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpperLower );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();

    /* Matrices */
    CHAMELEON_Complex64_t *A;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA * N );

    /* Set to alpha/beta */
    testing_start( &test_data );
    hres = CHAMELEON_zlaset( uplo, M, N, alpha, beta, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, 0.0 );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zset_std( args, uplo, alpha, beta, M, N, A, LDA );
    }

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zlaset;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zlaset_params[] = { "uplo", "m", "n", "lda", "alpha", "beta", NULL };
#else
const char *zlaset_params[] = { "mtxfmt", "nb", "uplo", "m", "n", "lda", "alpha", "beta", NULL };
#endif
const char *zlaset_output[] = { NULL };
const char *zlaset_outchk[] = { "||A||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlaset_init( void ) __attribute__( ( constructor ) );
void
testing_zlaset_init( void )
{
    test_zlaset.name   = "zlaset";
    test_zlaset.helper = "Initialize matrix with alpha (offdiag) and beta (diag)";
    test_zlaset.params = zlaset_params;
    test_zlaset.output = zlaset_output;
    test_zlaset.outchk = zlaset_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zlaset.fptr_desc = NULL;
#else
    test_zlaset.fptr_desc = testing_zlaset_desc;
#endif
    test_zlaset.fptr_std  = testing_zlaset_std;
    test_zlaset.next      = NULL;

    testing_register( &test_zlaset );
}
