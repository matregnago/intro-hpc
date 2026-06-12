/**
 *
 * @file testing_zpotri.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotri testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
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
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/lapacke.h>
#endif
#include <assert.h>

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zpotri_desc( run_arg_list_t *args, int check )
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

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Initialise the matrix with the random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    assert( hres == 0 );

    /* Calculates the inversed matrix */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zpotri_Tile_Async( uplo, descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zpotri_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotri( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the inverse */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamHermitian, uplo, ChamNonUnit, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zpotri_std( run_arg_list_t *args, int check )
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

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Initialise the matrix with the random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );

#if defined(CHAMELEON_TESTINGS_VENDOR)
    hres = LAPACKE_zpotrf( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
#else
    hres = CHAMELEON_zpotrf( uplo, N, A, LDA );
#endif
    assert( hres == 0 );

    /* Calculates the inversed matrix */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres += LAPACKE_zpotri( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotri( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres += CHAMELEON_zpotri( uplo, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zpotri( CblasColMajor, chameleon_lapack_const(uplo), N, A, LDA );
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
    testing_stop( &test_data, flops_zpotri( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the inverse */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_ztrtri_std( args, ChamHermitian, uplo, ChamNonUnit, N, A0, A, LDA );

        free( A0 );
    }
#endif

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zpotri;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zpotri_params[] = { "uplo", "n", "lda", "seedA", NULL };
#else
const char *zpotri_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
#endif
const char *zpotri_output[] = { NULL };
const char *zpotri_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotri_init( void ) __attribute__( ( constructor ) );
void
testing_zpotri_init( void )
{
    test_zpotri.name   = "zpotri";
    test_zpotri.helper = "Hermitian positive definite matrix inversion";
    test_zpotri.params = zpotri_params;
    test_zpotri.output = zpotri_output;
    test_zpotri.outchk = zpotri_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zpotri.fptr_desc = NULL;
#else
    test_zpotri.fptr_desc = testing_zpotri_desc;
#endif
    test_zpotri.fptr_std  = testing_zpotri_std;
    test_zpotri.next   = NULL;

    testing_register( &test_zpotri );
}
