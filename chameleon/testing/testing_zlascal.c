/**
 *
 * @file testing_zlascal.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
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

static cham_fixdbl_t
flops_zlascal( cham_uplo_t uplo, int _M, int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t minMN = (cham_fixdbl_t)chameleon_min( _M, _N );
    cham_fixdbl_t M     = _M;
    cham_fixdbl_t N     = _N;

    switch ( uplo ) {
        case ChamUpper:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + M * (cham_fixdbl_t)chameleon_max( 0, _N - _M );
            break;
        case ChamLower:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + N * (cham_fixdbl_t)chameleon_max( 0, _M - _N );
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* 1 multiplications per element */
    flops *= 6.;
#endif

    return flops;
}

int
testing_zlascal_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    int                   P     = parameters_getvalue_int( "P" );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", M );
    CHAMELEON_Complex64_t alpha = run_arg_get_complex64( args, "alpha", 1. );
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descAinit;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Scales the matrix */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlascal_Tile_Async( uplo, alpha, descA,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlascal_Tile( uplo, alpha, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlascal( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descAinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descAinit, seedA );

        hres += check_zscale( args, uplo, alpha, descAinit, descA );

        CHAMELEON_Desc_Destroy( &descAinit );
    }

    parameters_desc_destroy( &descA );

    return hres;
}

int
testing_zlascal_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   nb    = run_arg_get_nb( args );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", M );
    CHAMELEON_Complex64_t alpha = run_arg_get_complex64( args, "alpha", 1. );
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Scales the matrix */
    testing_start( &test_data );
    hres = CHAMELEON_zlascal( uplo, M, N, alpha, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlascal( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Ainit = malloc( sizeof(CHAMELEON_Complex64_t) * LDA * N );
        CHAMELEON_zplrnt( M, N, Ainit, LDA, seedA );

        hres += check_zscale_std( args, uplo, M, N, alpha, Ainit, A, LDA );

        free( Ainit );
    }

    free( A );

    return hres;
}

testing_t   test_zlascal;
const char *zlascal_params[] = { "mtxfmt", "nb", "uplo", "m", "n", "lda", "alpha", "seedA", NULL };
const char *zlascal_output[] = { NULL };
const char *zlascal_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlascal_init( void ) __attribute__( ( constructor ) );
void
testing_zlascal_init( void )
{
    test_zlascal.name   = "zlascal";
    test_zlascal.helper = "General matrix scaling";
    test_zlascal.params = zlascal_params;
    test_zlascal.output = zlascal_output;
    test_zlascal.outchk = zlascal_outchk;
    test_zlascal.fptr_desc = testing_zlascal_desc;
    test_zlascal.fptr_std  = testing_zlascal_std;
    test_zlascal.next   = NULL;

    testing_register( &test_zlascal );
}
