/**
 *
 * @file testing_zpotrs.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrs testing
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
#include <assert.h>

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zpotrs_desc( run_arg_list_t *args, int check )
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

    hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    assert( hres == 0 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zpotrs_Tile_Async( uplo, descA, descX,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zpotrs_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

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
testing_zpotrs_std( run_arg_list_t *args, int check )
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
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

#if defined(CHAMELEON_TESTINGS_VENDOR)
    hres = LAPACKE_zpotrf( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
#else
    hres = CHAMELEON_zpotrf( uplo, N, A, LDA );
#endif
    assert( hres == 0 );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres += LAPACKE_zpotrs( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, NRHS, A, LDA, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres += CHAMELEON_zpotrs( uplo, N, NRHS, A, LDA, X, LDB );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zpotrs( CblasColMajor, chameleon_lapack_const(uplo), N, NRHS, A, LDA, X, LDB );
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
    testing_stop( &test_data, flops_zpotrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        CHAMELEON_zplghe( (double)N, ChamUpperLower, N, A0, LDA, seedA );
        CHAMELEON_zplrnt( N, NRHS, B, LDB, seedB );

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

testing_t   test_zpotrs;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zpotrs_params[] = { "uplo",  "n",     "nrhs",
                                "lda",    "ldb", "seedA", "seedB", NULL };
#else
const char *zpotrs_params[] = { "mtxfmt", "nb",  "uplo",  "n",     "nrhs",
                                "lda",    "ldb", "seedA", "seedB", NULL };
#endif
const char *zpotrs_output[] = { NULL };
const char *zpotrs_outchk[] = { "||A||", "||A-fact(A)||", "||X||", "||B||", "||Ax-b||", "||Ax-b||/N/eps/(||A||||x||+||b||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrs_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrs_init( void )
{
    test_zpotrs.name   = "zpotrs";
    test_zpotrs.helper = "Hermitian positive definite solve (Cholesky)";
    test_zpotrs.params = zpotrs_params;
    test_zpotrs.output = zpotrs_output;
    test_zpotrs.outchk = zpotrs_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zpotrs.fptr_desc = NULL;
#else
    test_zpotrs.fptr_desc = testing_zpotrs_desc;
#endif
    test_zpotrs.fptr_std  = testing_zpotrs_std;
    test_zpotrs.next   = NULL;

    testing_register( &test_zpotrs );
}
