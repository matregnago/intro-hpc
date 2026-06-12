/**
 *
 * @file testing_zpotrf.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
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
#include <coreblas/lapacke.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zpotrf_desc( run_arg_list_t *args, int check )
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
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zpotrf_Tile_Async( uplo, descA,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, ChamUpperLower, descA0, seedA );

        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zpotrf_std( run_arg_list_t *args, int check )
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

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zpotrf( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zpotrf( uplo, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zpotrf( CblasColMajor, chameleon_lapack_const(uplo), N, A, LDA );
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
    testing_stop( &test_data, flops_zpotrf( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamHermitian, uplo, N, N, A0, A, LDA );

        free( A0 );
    }
#endif

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zpotrf;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zpotrf_params[] = { "uplo", "n", "lda", "seedA", NULL };
#else
const char *zpotrf_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
#endif
const char *zpotrf_output[] = { NULL };
const char *zpotrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrf_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrf_init( void )
{
    test_zpotrf.name   = "zpotrf";
    test_zpotrf.helper = "Hermitian positive definite factorization (Cholesky)";
    test_zpotrf.params = zpotrf_params;
    test_zpotrf.output = zpotrf_output;
    test_zpotrf.outchk = zpotrf_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zpotrf.fptr_desc = NULL;
#else
    test_zpotrf.fptr_desc = testing_zpotrf_desc;
#endif
    test_zpotrf.fptr_std  = testing_zpotrf_std;
    test_zpotrf.next   = NULL;

    testing_register( &test_zpotrf );
}
