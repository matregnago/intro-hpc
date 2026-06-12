/**
 *
 * @file testing_zsytrs.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrs testing
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
#include <assert.h>

int
testing_zsytrs_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );
    parameters_desc_create( "X", &descX, ChamComplexDouble, nb, nb, LDB, NRHS, N, NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    hres = CHAMELEON_zsytrf_Tile( uplo, descA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zsytrs_Tile_Async( uplo, descA, descX,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zsytrs_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, 0 /*flops_zsytrs( N, NRHS )*/ );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamSymmetric, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}

int
testing_zsytrs_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy( (double)N, uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    hres = CHAMELEON_zsytrf( uplo, N, A, LDA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    hres += CHAMELEON_zsytrs( uplo, N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, 0 /*flops_zsytrs( N, NRHS )*/ );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        CHAMELEON_zplgsy( (double)N, uplo, N, A0, LDA, seedA );
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        hres += check_zsolve_std( args, ChamSymmetric, ChamNoTrans, uplo, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    free( X );

    return hres;
}

testing_t   test_zsytrs;
const char *zsytrs_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                                "lda",    "ldb", "seedA", "seedB", NULL };
const char *zsytrs_output[] = { NULL };
const char *zsytrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsytrs_init( void ) __attribute__( ( constructor ) );
void
testing_zsytrs_init( void )
{
    test_zsytrs.name   = "zsytrs";
    test_zsytrs.helper = "Symmetric triangular solve";
    test_zsytrs.params = zsytrs_params;
    test_zsytrs.output = zsytrs_output;
    test_zsytrs.outchk = zsytrs_outchk;
    test_zsytrs.fptr_desc = testing_zsytrs_desc;
    test_zsytrs.fptr_std  = testing_zsytrs_std;
    test_zsytrs.next   = NULL;

    testing_register( &test_zsytrs );
}
