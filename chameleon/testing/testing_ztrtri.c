/**
 *
 * @file testing_ztrtri.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri testing
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
#endif

int
testing_ztrtri_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the inversed matrices */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_ztrtri_Tile_Async( uplo, diag, descA,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_ztrtri_Tile( uplo, diag, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrtri( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the inverse */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamTriangular, uplo, diag, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

int
testing_ztrtri_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         api   = parameters_getvalue_int( "api" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );

    /* Calculates the inversed matrices */
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_ztrtri( uplo, diag, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_ztrtri( CblasColMajor, chameleon_lapack_const(uplo), chameleon_lapack_const(diag), N, A, LDA );
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
    testing_stop( &test_data, flops_ztrtri( N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the inverse */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_ztrtri_std( args, ChamTriangular, uplo, diag, N, A0, A, LDA );

        free( A0 );
    }

    free( A );

    return hres;
}

testing_t   test_ztrtri;
const char *ztrtri_params[] = { "mtxfmt", "nb", "uplo", "diag", "n", "lda", "seedA", NULL };
const char *ztrtri_output[] = { NULL };
const char *ztrtri_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrtri_init( void ) __attribute__( ( constructor ) );
void
testing_ztrtri_init( void )
{
    test_ztrtri.name   = "ztrtri";
    test_ztrtri.helper = "Triangular matrix inversion";
    test_ztrtri.params = ztrtri_params;
    test_ztrtri.output = ztrtri_output;
    test_ztrtri.outchk = ztrtri_outchk;
    test_ztrtri.fptr_desc = testing_ztrtri_desc;
    test_ztrtri.fptr_std  = testing_ztrtri_std;
    test_ztrtri.next   = NULL;

    testing_register( &test_ztrtri );
}
