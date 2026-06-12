/**
 *
 * @file testing_zunglq.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunglq testing
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
testing_zunglq_desc( run_arg_list_t *args, int check )
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
    int      K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      RH     = run_arg_get_int( args, "qra", 0 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descT, *descQ;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( M > N ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for unglq (M > N)\n" );
        }
        return -1;
    }

    if ( K > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for unglq (K > M)\n" );
        }
        return -1;
    }

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, K, N, P, Q );
    CHAMELEON_Desc_Create(
        &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    hres = CHAMELEON_zgelqf_Tile( descA, descT );
    if ( hres != CHAMELEON_SUCCESS ) {
        CHAMELEON_Desc_Destroy( &descA );
        CHAMELEON_Desc_Destroy( &descT );
        CHAMELEON_Desc_Destroy( &descQ );
        return hres;
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zunglq_Tile_Async( descA, descT, descQ,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
        CHAMELEON_Desc_Flush( descQ, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zunglq_Tile( descA, descT, descQ );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunglq( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        hres += check_zortho( args, descQ );
        hres += check_zgelqf( args, descA0, descA, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descQ );

    return hres;
}

int
testing_zunglq_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int nb    = run_arg_get_nb( args );
    int ib    = run_arg_get_ib( args );
    int P     = parameters_getvalue_int( "P" );
    int N     = run_arg_get_int( args, "N", 1000 );
    int M     = run_arg_get_int( args, "M", N );
    int K     = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int LDA   = run_arg_get_int( args, "LDA", M );
    int RH    = run_arg_get_int( args, "qra", 0 );
    int seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *Qlap;
    CHAM_desc_t           *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( M > N ) {
        fprintf( stderr, "SKIPPED: Incorrect parameters for unglq (M > N)\n" );
        return -1;
    }

    if ( K > M ) {
        fprintf( stderr, "SKIPPED: Incorrect parameters for unglq (K > M)\n" );
        return -1;
    }

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    A    = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    Qlap = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    CHAMELEON_Alloc_Workspace_zgels( K, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( K, N, A, LDA, seedA );
    hres = CHAMELEON_zgelqf( K, N, A, LDA, descT );
    if ( hres != CHAMELEON_SUCCESS ) {
        free( A    );
        free( Qlap );
        CHAMELEON_Desc_Destroy( &descT );
        return hres;
    }

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zunglq( M, N, K, A, LDA, descT, Qlap, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunglq( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zplrnt( K, N, A0, LDA, seedA );

        hres += check_zortho_std( args, M, N, Qlap, LDA );
        hres += check_zgelqf_std( args, M, N, K, A0, A, LDA, Qlap, LDA );

        free( A0 );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descT );
    free( Qlap );

    return hres;
}

testing_t   test_zunglq;
const char *zunglq_params[] = { "mtxfmt", "nb", "ib", "m", "n", "k", "lda", "qra", "seedA", NULL };
const char *zunglq_output[] = { NULL };
const char *zunglq_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunglq_init( void ) __attribute__( ( constructor ) );
void
testing_zunglq_init( void )
{
    test_zunglq.name   = "zunglq";
    test_zunglq.helper = "Q generation (LQ)";
    test_zunglq.params = zunglq_params;
    test_zunglq.output = zunglq_output;
    test_zunglq.outchk = zunglq_outchk;
    test_zunglq.fptr_desc = testing_zunglq_desc;
    test_zunglq.fptr_std  = testing_zunglq_std;
    test_zunglq.next   = NULL;

    testing_register( &test_zunglq );
}
