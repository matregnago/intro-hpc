/**
 *
 * @file testing_zlauum.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
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
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#endif

static cham_fixdbl_t
flops_zlauum( int N )
{
    cham_fixdbl_t flops = flops_zpotri( N ) - flops_ztrtri( N );
    return flops;
}

int
testing_zlauum_desc( run_arg_list_t *args, int check )
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

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe_Tile( 0., uplo, descA, seedA );

    /* Calculates the matrix product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlauum_Tile_Async( uplo, descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlauum_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlauum( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( 0., uplo, descA0, seedA );

        hres += check_zlauum( args, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

int
testing_zlauum_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int         api   = parameters_getvalue_int( "api" );
#endif
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

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe( 0., uplo, N, A, LDA, seedA );

    /* Calculates the matrix product */
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zlauum( uplo, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zlauum( CblasColMajor, chameleon_lapack_const(uplo), N, A, LDA );
        break;
#endif
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlauum( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplghe( 0., uplo, N, A0, LDA, seedA );

        hres += check_zlauum_std( args, uplo, N, A0, A, LDA );

        free( A0 );
    }

    free( A );

    return hres;
}

testing_t   test_zlauum;
const char *zlauum_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
const char *zlauum_output[] = { NULL };
const char *zlauum_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlauum_init( void ) __attribute__( ( constructor ) );
void
testing_zlauum_init( void )
{
    test_zlauum.name   = "zlauum";
    test_zlauum.helper = "Triangular in-place matrix-matrix computation for Cholesky inversion";
    test_zlauum.params = zlauum_params;
    test_zlauum.output = zlauum_output;
    test_zlauum.outchk = zlauum_outchk;
    test_zlauum.fptr_desc = testing_zlauum_desc;
    test_zlauum.fptr_std  = testing_zlauum_std;
    test_zlauum.next   = NULL;

    testing_register( &test_zlauum );
}
