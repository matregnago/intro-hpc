/**
 *
 * @file testing_zgepdf_qdwh.c
 *
 * @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2020-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qdwh testing
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
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
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#endif

int
testing_zgepdf_qdwh_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      LDB    = run_arg_get_int( args, "LDB", N );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );
    double   cond   = run_arg_get_double( args, "cond", 1.e16 );
    int      mode   = run_arg_get_int( args, "mode", 4 );
    int      runtime;

    /* Descriptors */
    CHAM_desc_t *descA, *descH, *descA0;
    gepdf_info_t info = { .itQR = 0, .itPO = 0, .flops = 0. };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    CHAMELEON_Get( CHAMELEON_RUNTIME, &runtime );
    if ( runtime == RUNTIME_SCHED_PARSEC ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The QDWH polar decomposition is not supported with PaRSEC\n" );
        }
        return -1;
    }

    if ( N > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The QDWH polar decomposition is performed only when M >= N\n" );
        }
        return -1;
    }

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descH, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms_Tile(
        ChamDistUniform, seedA, ChamNonsymPosv, NULL, mode, cond, 1., descA );
    if ( hres != 0 ) {
        return hres;
    }
    /*
     * descA0 is defined here because of the cost of zlatms. To copy descA in descA0
     * now prevents to call it again later in the check (indeed descA is modified
     * with the call to CHAMELEON_zgepdf_qdwh_Tile[_Async]).
     */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zlacpy_Tile( ChamUpperLower, descA, descA0 );
    }

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgepdf_qdwh_Tile_Async( descA, descH, &info,
                                                 test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descH, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgepdf_qdwh_Tile( descA, descH, &info );
    }
    test_data.hres = hres;
    testing_stop( &test_data, info.flops );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zxxpd( args, descA0, descA, descH );
        hres += check_zortho( args, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descH );

    return hres;
}

int
testing_zgepdf_qdwh_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      LDB    = run_arg_get_int( args, "LDB", N );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    double   cond   = run_arg_get_double( args, "cond", 1.e16 );
    int      mode   = run_arg_get_int( args, "mode", 4 );
    int      runtime;

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *H, *A0;
    gepdf_info_t info = { .itQR = 0, .itPO = 0, .flops = 0. };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    CHAMELEON_Get( CHAMELEON_RUNTIME, &runtime );
    if ( runtime == RUNTIME_SCHED_PARSEC ) {
        fprintf( stderr, "SKIPPED: The QDWH polar decomposition is not supported with PaRSEC\n" );
        return -1;
    }

    if ( N > M ) {
        fprintf( stderr, "SKIPPED: The QDWH polar decomposition is performed only when M >= N\n" );
        return -1;
    }

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    H = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms( M, N, ChamDistUniform, seedA, ChamNonsymPosv, NULL, mode, cond, 1., A, LDA );
    if ( hres != 0 ) {
        free( A );
        free( H );
        return hres;
    }
    /*
     * A0 is defined here because of the cost of zlatms. To copy A in A0
     * now prevents to call it again later in the check (indeed A is modified
     * with the call to CHAMELEON_zgepdf_qdwh).
     */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zlacpy( ChamUpperLower, M, N, A, LDA, A0, LDA );
    }

    /* Calculates the norm */
    testing_start( &test_data );
    hres = CHAMELEON_zgepdf_qdwh( M, N, A, LDA, H, LDB, &info );
    test_data.hres = hres;
    testing_stop( &test_data, info.flops );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zxxpd_std ( args, M, N, A0, A, LDA, H, LDB );
        hres += check_zortho_std( args, M, N, A, LDA );

        free( A0 );
    }

    free( A );
    free( H );

    return hres;
}

testing_t   test_zgepdf_qdwh;
const char *zgepdf_qdwh_params[] = { "mtxfmt", "nb",    "ib",   "m",    "n", "lda",
                                     "ldb",    "seedA", "cond", "mode", NULL };
const char *zgepdf_qdwh_output[] = { NULL };
const char *zgepdf_qdwh_outchk[] = { "||A||", "||A-fact(A)||", "||I-QQ'||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgepdf_qdwh_init( void ) __attribute__( ( constructor ) );
void
testing_zgepdf_qdwh_init( void )
{
    test_zgepdf_qdwh.name   = "zgepdf_qdwh";
    test_zgepdf_qdwh.helper = "Polar decomposition factorization with QDWH algorithm";
    test_zgepdf_qdwh.params = zgepdf_qdwh_params;
    test_zgepdf_qdwh.output = zgepdf_qdwh_output;
    test_zgepdf_qdwh.outchk = zgepdf_qdwh_outchk;
    test_zgepdf_qdwh.fptr_desc = testing_zgepdf_qdwh_desc;
    test_zgepdf_qdwh.fptr_std  = testing_zgepdf_qdwh_std;
    test_zgepdf_qdwh.next   = NULL;

    testing_register( &test_zgepdf_qdwh );
}
