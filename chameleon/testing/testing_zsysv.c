/**
 *
 * @file testing_zsysv.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsysv testing
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

static cham_fixdbl_t
flops_zsysv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

int
testing_zsysv_desc( run_arg_list_t *args, int check )
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

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsysv_Tile_Async( uplo, descA, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsysv_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsysv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplgsy_Tile( (double)N, ChamUpperLower, descA0, seedA );

        hres += check_zxxtrf( args, ChamSymmetric, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamSymmetric, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}

int
testing_zsysv_std( run_arg_list_t *args, int check )
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

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zsysv( uplo, N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsysv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0, *B;

        /* Check the factorization */
        A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplgsy( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamSymmetric, uplo, N, N, A0, A, LDA );

        /* Check the solve */
        B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        CHAMELEON_zplgsy( (double)N, uplo, N, A0, LDA, seedA );
        hres += check_zsolve_std( args, ChamSymmetric, ChamNoTrans, uplo, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    free( X );

    return hres;
}

testing_t   test_zsysv;
const char *zsysv_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                               "lda",    "ldb", "seedA", "seedB", NULL };
const char *zsysv_output[] = { NULL };
const char *zsysv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsysv_init( void ) __attribute__( ( constructor ) );
void
testing_zsysv_init( void )
{
    test_zsysv.name   = "zsysv";
    test_zsysv.helper = "Symmetrix linear system solve";
    test_zsysv.params = zsysv_params;
    test_zsysv.output = zsysv_output;
    test_zsysv.outchk = zsysv_outchk;
    test_zsysv.fptr_desc = testing_zsysv_desc;
    test_zsysv.fptr_std  = testing_zsysv_std;
    test_zsysv.next   = NULL;

    testing_register( &test_zsysv );
}
