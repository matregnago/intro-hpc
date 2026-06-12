/**
 *
 * @file testing_zgetrs_nopiv.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrs_nopiv testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
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
#include <assert.h>

int
testing_zgetrs_nopiv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int                   LDA   = run_arg_get_int( args, "LDA", N );
    int                   LDB   = run_arg_get_int( args, "LDB", N );
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    CHAMELEON_Complex64_t bump  = run_arg_get_complex64( args, "bump", (CHAMELEON_Complex64_t)N );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );
    parameters_desc_create( "X", &descX, ChamComplexDouble, nb, nb, LDB, NRHS, N, NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA, seedA   );
    CHAMELEON_zplgtr_Tile( bump, ChamLower, descA, seedA+1 );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    hres = CHAMELEON_zgetrf_nopiv_Tile( descA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zgetrs_nopiv_Tile_Async( descA, descX,
                                                   test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zgetrs_nopiv_Tile( descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA0, seedA   );
        CHAMELEON_zplgtr_Tile( bump, ChamLower, descA0, seedA+1 );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower,
                              descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}

int
testing_zgetrs_nopiv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   nb    = run_arg_get_nb( args );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int                   LDA   = run_arg_get_int( args, "LDA", N );
    int                   LDB   = run_arg_get_int( args, "LDB", N );
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    CHAMELEON_Complex64_t bump  = run_arg_get_complex64( args, "bump", (CHAMELEON_Complex64_t)N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr( 0,    ChamUpper, N, N, A, LDA, seedA   );
    CHAMELEON_zplgtr( bump, ChamLower, N, N, A, LDA, seedA+1 );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    hres = CHAMELEON_zgetrf_nopiv( N, N, A, LDA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    hres += CHAMELEON_zgetrs_nopiv( ChamNoTrans, N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        CHAMELEON_zplgtr( 0,    ChamUpper, N, N, A0, LDA, seedA   );
        CHAMELEON_zplgtr( bump, ChamLower, N, N, A0, LDA, seedA+1 );
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        hres += check_zsolve_std( args, ChamGeneral, ChamNoTrans, ChamUpperLower,
                                  N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }

    free( A );
    free( X );

    return hres;
}

testing_t   test_zgetrs_nopiv;
const char *zgetrs_nopiv_params[] = { "mtxfmt", "nb", "n", "nrhs", "lda", "ldb", "seedA", "seedB", "bump", NULL };
const char *zgetrs_nopiv_output[] = { NULL };
const char *zgetrs_nopiv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrs_nopiv_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrs_nopiv_init( void )
{
    test_zgetrs_nopiv.name   = "zgetrs_nopiv";
    test_zgetrs_nopiv.helper = "General triangular solve (LU without pivoting)";
    test_zgetrs_nopiv.params = zgetrs_nopiv_params;
    test_zgetrs_nopiv.output = zgetrs_nopiv_output;
    test_zgetrs_nopiv.outchk = zgetrs_nopiv_outchk;
    test_zgetrs_nopiv.fptr_desc = testing_zgetrs_nopiv_desc;
    test_zgetrs_nopiv.fptr_std  = testing_zgetrs_nopiv_std;
    test_zgetrs_nopiv.next   = NULL;

    testing_register( &test_zgetrs_nopiv );
}
