/**
 *
 * @file testing_zposv.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zposv testing
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

static cham_fixdbl_t
flops_zposv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zposv_desc( run_arg_list_t *args, int check )
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
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zposv_Tile_Async( uplo, descA, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zposv_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zposv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and the residual */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, ChamUpperLower, descA0, seedA );

        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamHermitian, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}
#endif

int
testing_zposv_std( run_arg_list_t *args, int check )
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
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N    );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zposv( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zposv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zposv( uplo, N, NRHS, A, LDA, X, LDB );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zposv( CblasColMajor, chameleon_lapack_const(uplo), N, NRHS, A, LDA, X, LDB );
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
    testing_stop( &test_data, flops_zposv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residual */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        /* Check the factorization */
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamHermitian, uplo, N, N, A0, A, LDA );

        /* Check the solve */
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );
        hres += check_zsolve_std( args, ChamHermitian, ChamNoTrans, uplo, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }
#endif

    free( A );
    free( X );

    (void)check;
    return hres;
}

testing_t   test_zposv;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zposv_params[] = { "uplo",  "n",     "nrhs",
                               "lda",    "ldb", "seedA", "seedB", NULL };
#else
const char *zposv_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                               "lda",    "ldb", "seedA", "seedB", NULL };
#endif
const char *zposv_output[] = { NULL };
const char *zposv_outchk[] = { "||A||", "||A-fact(A)||", "||X||", "||B||", "||Ax-b||", "||Ax-b||/N/eps/(||A||||x||+||b||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zposv_init( void ) __attribute__( ( constructor ) );
void
testing_zposv_init( void )
{
    test_zposv.name   = "zposv";
    test_zposv.helper = "Hermitian positive definite linear system solve (Cholesky)";
    test_zposv.params = zposv_params;
    test_zposv.output = zposv_output;
    test_zposv.outchk = zposv_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zposv.fptr_desc = NULL;
#else
    test_zposv.fptr_desc = testing_zposv_desc;
#endif
    test_zposv.fptr_std  = testing_zposv_std;
    test_zposv.next   = NULL;

    testing_register( &test_zposv );
}
