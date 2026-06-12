/**
 *
 * @file testing_zsytrf.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zsytrf_desc( run_arg_list_t *args, int check )
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

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsytrf_Tile_Async( uplo, descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsytrf_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplgsy_Tile( (double)N, ChamUpperLower, descA0, seedA );

        hres += check_zxxtrf( args, ChamSymmetric, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

int
testing_zsytrf_std( run_arg_list_t *args, int check )
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

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy( (double)N, uplo, N, A, LDA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zsytrf( uplo, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplgsy( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamSymmetric, uplo, N, N, A0, A, LDA );

        free( A0 );
    }

    free( A );

    return hres;
}

testing_t   test_zsytrf;
const char *zsytrf_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
const char *zsytrf_output[] = { NULL };
const char *zsytrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsytrf_init( void ) __attribute__( ( constructor ) );
void
testing_zsytrf_init( void )
{
    test_zsytrf.name   = "zsytrf";
    test_zsytrf.helper = "Symmetric trinagular factorization";
    test_zsytrf.params = zsytrf_params;
    test_zsytrf.output = zsytrf_output;
    test_zsytrf.outchk = zsytrf_outchk;
    test_zsytrf.fptr_desc = testing_zsytrf_desc;
    test_zsytrf.fptr_std  = testing_zsytrf_std;
    test_zsytrf.next   = NULL;

    testing_register( &test_zsytrf );
}
