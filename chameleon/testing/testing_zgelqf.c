/**
 *
 * @file testing_zgelqf.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
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

int
testing_zgelqf_desc( run_arg_list_t *args, int check )
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
    int      RH     = run_arg_get_int( args, "qra", 4 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgelqf_Tile_Async( descA, descT, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgelqf_Tile( descA, descT );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgelqf( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descQ;
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, N, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zunglq_Tile( descA, descT, descQ );

        hres += check_zgelqf( args, descA0, descA, descQ );
        hres += check_zortho( args, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

int
testing_zgelqf_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int nb    = run_arg_get_nb( args );
    int ib    = run_arg_get_ib( args );
    int P     = parameters_getvalue_int( "P" );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int LDA   = run_arg_get_int( args, "LDA", M );
    int RH    = run_arg_get_int( args, "qra", 4 );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int Q     = parameters_compute_q( P );
    int K     = chameleon_min( M, N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;
    CHAM_desc_t           *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgelqf( M, N, A, LDA, descT );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgelqf( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Qlap = malloc( sizeof(CHAMELEON_Complex64_t) * N*N );
        CHAMELEON_Complex64_t *A0   = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

        CHAMELEON_zplrnt( M, N, A0, LDA, seedA );

        CHAMELEON_zunglq( N, N, K, A, LDA, descT, Qlap, N );

        hres += check_zgelqf_std( args, N, N, K, A0, A, LDA, Qlap, N );
        hres += check_zortho_std( args, N, N, Qlap, N );

        free( A0 );
        free( Qlap );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

testing_t   test_zgelqf;
const char *zgelqf_params[] = { "mtxfmt", "nb", "ib", "m", "n", "lda", "qra", "seedA", NULL };
const char *zgelqf_output[] = { NULL };
const char *zgelqf_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgelqf_init( void ) __attribute__( ( constructor ) );
void
testing_zgelqf_init( void )
{
    test_zgelqf.name   = "zgelqf";
    test_zgelqf.helper = "General LQ factorization";
    test_zgelqf.params = zgelqf_params;
    test_zgelqf.output = zgelqf_output;
    test_zgelqf.outchk = zgelqf_outchk;
    test_zgelqf.fptr_desc = testing_zgelqf_desc;
    test_zgelqf.fptr_std  = testing_zgelqf_std;
    test_zgelqf.next   = NULL;

    testing_register( &test_zgelqf );
}
