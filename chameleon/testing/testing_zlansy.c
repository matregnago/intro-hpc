/**
 *
 * @file testing_zlansy.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Philippe Swartvagher
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

static cham_fixdbl_t
flops_zlansy( cham_normtype_t ntype, int _N )
{
    cham_fixdbl_t flops   = 0.;
    cham_fixdbl_t coefabs = 1.;
    cham_fixdbl_t N       = _N;
    cham_fixdbl_t size    = ( N * ( N + 1 ) ) / 2.;
#if defined(PRECISION_z) || defined(PRECISION_c)
    coefabs = 3.;
#endif

    switch ( ntype ) {
        case ChamMaxNorm:
            flops = coefabs * size;
            break;
        case ChamOneNorm:
        case ChamInfNorm:
            flops = coefabs * size + N * ( N - 1 );
            break;
        case ChamFrobeniusNorm:
            flops = ( coefabs + 1. ) * size;
            break;
        default:;
    }
    (void)flops;
    return sizeof( CHAMELEON_Complex64_t ) * size;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zlansy_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async     = parameters_getvalue_int( "async" );
    int                   nb        = run_arg_get_nb( args );
    cham_normtype_t       norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t           uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N         = run_arg_get_int( args, "N", 1000 );
    int                   LDA       = run_arg_get_int( args, "LDA", N );
    int                   seedA     = run_arg_get_int( args, "seedA", testing_ialea() );
    CHAMELEON_Complex64_t bump      = testing_zalea();

    /* Descriptors */
    double       norm;
    CHAM_desc_t *descA;

    bump = run_arg_get_complex64( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( bump, uplo, descA, seedA );

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlansy_Tile_Async( norm_type, uplo, descA, &norm,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        norm = CHAMELEON_zlansy_Tile( norm_type, uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlansy( norm_type, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_znorm( args, ChamSymmetric, norm_type, uplo, ChamNonUnit, norm, descA );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zlansy_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int                   api       = parameters_getvalue_int( "api" );
#endif
    int                   nb        = run_arg_get_nb( args );
    cham_normtype_t       norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t           uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N         = run_arg_get_int( args, "N", 1000 );
    int                   LDA       = run_arg_get_int( args, "LDA", N );
    int                   seedA     = run_arg_get_int( args, "seedA", testing_ialea() );
    CHAMELEON_Complex64_t bump      = testing_zalea();

    /* Descriptors */
    double norm;
    CHAMELEON_Complex64_t *A;

    bump = run_arg_get_complex64( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy( bump, uplo, N, A, LDA, seedA );

    /* Calculates the norm */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    norm = LAPACKE_zlansy( LAPACK_COL_MAJOR, chameleon_lapack_const( norm_type ), uplo, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlansy( norm_type, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        norm = CHAMELEON_zlansy( norm_type, uplo, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        norm = CHAMELEON_lapacke_zlansy( CblasColMajor, chameleon_lapack_const( norm_type ), chameleon_lapack_const(uplo), N, A, LDA );
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
    testing_stop( &test_data, flops_zlansy( norm_type, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_znorm_std( args, ChamSymmetric, norm_type, uplo, ChamNonUnit, norm, N, N, A, LDA );
    }
#endif

    free( A );

    (void)norm;
    (void)check;
    return hres;
}

testing_t   test_zlansy;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zlansy_params[] = { "norm", "uplo", "n", "lda", "seedA", "bump", NULL };
#else
const char *zlansy_params[] = { "mtxfmt", "nb", "norm", "uplo", "n", "lda", "seedA", "bump", NULL };
#endif
const char *zlansy_output[] = { NULL };
const char *zlansy_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlansy_init( void ) __attribute__( ( constructor ) );
void
testing_zlansy_init( void )
{
    test_zlansy.name   = "zlansy";
    test_zlansy.helper = "Symmetric matrix norm";
    test_zlansy.params = zlansy_params;
    test_zlansy.output = zlansy_output;
    test_zlansy.outchk = zlansy_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zlansy.fptr_desc = NULL;
#else
    test_zlansy.fptr_desc = testing_zlansy_desc;
#endif
    test_zlansy.fptr_std  = testing_zlansy_std;
    test_zlansy.next   = NULL;

    testing_register( &test_zlansy );
}
