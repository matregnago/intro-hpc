/**
 *
 * @file testing_zlacpy.c
 *
 * @copyright 2019-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy testing
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
flops_zlacpy( cham_uplo_t uplo, int _M, int _N )
{
    cham_fixdbl_t flops;
    cham_fixdbl_t M = _M;
    cham_fixdbl_t N = _N;

    switch ( uplo ) {
        case ChamUpper:
            if ( N > M ) {
                flops = ( M * ( M + 1 ) / 2 ) + M * ( N - M );
            }
            else {
                flops = N * ( N + 1 ) / 2;
            }
            break;
        case ChamLower:
            if ( M > N ) {
                flops = ( N * ( N + 1 ) / 2 ) + N * ( M - N );
            }
            else {
                flops = M * ( M + 1 ) / 2;
            }
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }
    flops *= sizeof( CHAMELEON_Complex64_t );

    return flops;
}

int
testing_zlacpy_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         LDB   = run_arg_get_int( args, "LDB", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA, *descB;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, N, M, N );

    /* Fills each matrix with different random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    /* We use seedA + 1, just to create a variation in B */
    CHAMELEON_zplrnt_Tile( descB, seedA + 1 );

    /* Makes a copy of descA to descB */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlacpy_Tile_Async( uplo, descA, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlacpy_Tile( uplo, descA, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlacpy( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks their differences */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zmatrices( args, uplo, descA, descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}

int
testing_zlacpy_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         api   = parameters_getvalue_int( "api" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         LDB   = run_arg_get_int( args, "LDB", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *B;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N );

    /* Fills each matrix with different random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );
    /* We use seedA + 1, just to create a variation in B */
    CHAMELEON_zplrnt( M, N, B, LDB, seedA + 1 );

    /* Makes a copy of descA to descB */
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zlacpy( uplo, M, N, A, LDA, B, LDB );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zlacpy( CblasColMajor, chameleon_lapack_const(uplo), M, N, A, LDA, B, LDB );
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
    testing_stop( &test_data, flops_zlacpy( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks their differences */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zmatrices_std( args, uplo, M, N, A, LDA, B, LDB );
    }

    free( A );
    free( B );

    return hres;
}

testing_t   test_zlacpy;
const char *zlacpy_params[] = { "mtxfmt", "nb", "uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zlacpy_output[] = { NULL };
const char *zlacpy_outchk[] = { "||A||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlacpy_init( void ) __attribute__( ( constructor ) );
void
testing_zlacpy_init( void )
{
    test_zlacpy.name   = "zlacpy";
    test_zlacpy.helper = "General matrix copy";
    test_zlacpy.params = zlacpy_params;
    test_zlacpy.output = zlacpy_output;
    test_zlacpy.outchk = zlacpy_outchk;
    test_zlacpy.fptr_desc = testing_zlacpy_desc;
    test_zlacpy.fptr_std  = testing_zlacpy_std;
    test_zlacpy.next   = NULL;

    testing_register( &test_zlacpy );
}
