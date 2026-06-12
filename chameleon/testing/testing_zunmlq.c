/**
 *
 * @file testing_zunmlq.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq testing
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
#include <assert.h>

int
testing_zunmlq_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_nb( args );
    int          ib     = run_arg_get_ib( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_side_t  side   = run_arg_get_side( args, "side", ChamLeft );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int          LDA    = run_arg_get_int( args, "LDA", K );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    int          RH     = run_arg_get_int( args, "qra", 4 );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedC  = run_arg_get_int( args, "seedC", testing_ialea() );
    int          Q      = parameters_compute_q( P );

    /* Descriptors */
    int          An;
    CHAM_desc_t *descA, *descT, *descC;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Calculates the dimensions according to the transposition and the side */
    An = ( side == ChamLeft ) ? M : N;

    if ( K > An )
    {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameter K for unmqr\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, An, 0, 0, K, An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( K, An, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgelqf_Tile( descA, descT );
    assert( hres == 0 );

    /* Computes unmlq */
    testing_start( &test_data );
    if ( async ) {
        hres += CHAMELEON_zunmlq_Tile_Async( side, trans, descA, descT, descC,
                                             test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres += CHAMELEON_zunmlq_Tile( side, trans, descA, descT, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunmlq( side, M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descC0 = CHAMELEON_Desc_Copy( descC, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descQ;

        CHAMELEON_zplrnt_Tile( descC0, seedC );

        CHAMELEON_Desc_Create(
            &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, An, An, 0, 0, An, An, P, Q );
        CHAMELEON_zunglq_Tile( descA, descT, descQ );

        hres += check_zqc( args, side, trans, descC0, descQ, descC );

        CHAMELEON_Desc_Destroy( &descC0 );
        CHAMELEON_Desc_Destroy( &descQ );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

int
testing_zunmlq_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int          nb    = run_arg_get_nb( args );
    int          ib    = run_arg_get_ib( args );
    int          P     = parameters_getvalue_int( "P" );
    cham_side_t  side  = run_arg_get_side( args, "side", ChamLeft );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          M     = run_arg_get_int( args, "M", N );
    int          K     = run_arg_get_int( args, "K", N );
    int          LDA   = run_arg_get_int( args, "LDA", K );
    int          LDC   = run_arg_get_int( args, "LDC", M );
    int          RH    = run_arg_get_int( args, "qra", 4 );
    int          seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    int          Q     = parameters_compute_q( P );

    /* Descriptors */
    int                    An;
    CHAMELEON_Complex64_t *A, *C;
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

    /* Calculates the dimensions according to the transposition and the side */
    An = ( side == ChamLeft ) ? M : N;

    if ( K > An )
    {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameter K for unmqr\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    C = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N  );
    CHAMELEON_Alloc_Workspace_zgels( K, An, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( K, An, A, LDA, seedA );
    CHAMELEON_zplrnt( M, N , C, LDC, seedC );

    /* Computes the factorization */
    hres = CHAMELEON_zgelqf( K, An, A, LDA, descT );

    /* Computes unmlq */
    testing_start( &test_data );
    hres += CHAMELEON_zunmlq( side, trans, M, N, K, A, LDA, descT, C, LDC );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zunmlq( side, M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and orthogonality */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *C0   = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );
        CHAMELEON_Complex64_t *Qlap = malloc( sizeof(CHAMELEON_Complex64_t) * An*An );

        CHAMELEON_zplrnt( M, N, C0, LDC, seedC );
        CHAMELEON_zunglq( An, An, K, A, LDA, descT, Qlap, An );

        hres += check_zqc_std( args, side, trans, M, N, C0, C, LDC, Qlap, An );

        free( C0 );
        free( Qlap );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descT );
    free( C );

    return hres;
}

testing_t   test_zunmlq;
const char *zunmlq_params[] = { "mtxfmt", "nb",  "ib",  "side", "trans", "m",     "n",
                                "k",      "lda", "ldc", "qra",  "seedA", "seedC", NULL };
const char *zunmlq_output[] = { NULL };
const char *zunmlq_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zunmlq_init( void ) __attribute__( ( constructor ) );
void
testing_zunmlq_init( void )
{
    test_zunmlq.name   = "zunmlq";
    test_zunmlq.helper = "Q application (LQ)";
    test_zunmlq.params = zunmlq_params;
    test_zunmlq.output = zunmlq_output;
    test_zunmlq.outchk = zunmlq_outchk;
    test_zunmlq.fptr_desc = testing_zunmlq_desc;
    test_zunmlq.fptr_std  = testing_zunmlq_std;
    test_zunmlq.next   = NULL;

    testing_register( &test_zunmlq );
}
