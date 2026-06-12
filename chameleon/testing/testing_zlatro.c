/**
 *
 * @file testing_zlatro.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro testing
 *
 * @version 1.3.0
 * @author Brieuc Nicolas
 * @date 2025-12-17
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"

int
testing_zlatro_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                    async = parameters_getvalue_int( "async" );
    int                    nb    = run_arg_get_nb( args );
    cham_trans_t           trans = run_arg_get_trans( args, "trans", ChamTrans );
    cham_uplo_t            uplo  = run_arg_get_uplo( args, "uplo", ChamUpperLower );
    int                    N     = run_arg_get_int( args, "N", 1000 );
    int                    M     = run_arg_get_int( args, "M", N );
    int                    LDA   = run_arg_get_int( args, "LDA", M );
    int                    LDB   = run_arg_get_int( args, "LDB", N );
    unsigned long long int seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA, *descB;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    if( trans == ChamNoTrans ) {
        fprintf(stderr, "testing_zlatro_desc: Nothing to do for ChamNoTrans\n" );
        return CHAMELEON_SUCCESS;
    }

    /* Creates the matrix */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, M, N, M );

    /* Fills the matrix with random values */
    if( uplo == ChamUpperLower ) {
        CHAMELEON_zplrnt_Tile( descA, seedA );
    }
    else {
        CHAMELEON_zplgtr_Tile( 0., uplo, descA, seedA );
    }
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0, 0, descB );

    /* Randomly generates the matrix */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlatro_Tile_Async( uplo, trans, descA, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlatro_Tile( uplo, trans, descA, descB );
    }
    testing_stop( &test_data, 0 );
    test_data.hres = hres;

    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_ztranspose( args, uplo, trans, descA, descB );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}

int
testing_zlatro_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                    nb    = run_arg_get_nb( args );
    cham_trans_t           trans = run_arg_get_trans( args, "trans", ChamTrans );
    cham_uplo_t            uplo  = run_arg_get_uplo( args, "uplo", ChamUpperLower );
    int                    N     = run_arg_get_int( args, "N", 1000 );
    int                    M     = run_arg_get_int( args, "M", N );
    int                    LDA   = run_arg_get_int( args, "LDA", M );
    int                    LDB   = run_arg_get_int( args, "LDB", N );
    unsigned long long int seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *B;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    if( trans == ChamNoTrans ) {
        fprintf(stderr, "testing_zlatro_desc: Nothing to do for ChamNoTrans\n" );
        return CHAMELEON_SUCCESS;
    }

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA * N );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB * M );

    /* Fills the matrix with random values */
    if( uplo == ChamUpperLower ) {
        CHAMELEON_zplrnt( M, N, A, LDA, seedA );
    }
    else {
        CHAMELEON_zplgtr( 0., uplo, M, N, A, LDA, seedA );
    }
    CHAMELEON_zlaset( ChamUpperLower, N, M, 0., 0., B, LDB );

    /* Randomly generates the matrix */
    testing_start( &test_data );
    hres = CHAMELEON_zlatro( uplo, trans, M, N, A, LDA, B, LDB );
    testing_stop( &test_data, 0 );
    test_data.hres = hres;

    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_ztranspose_std( args, uplo, trans, M, N,
                                      A, LDA, B, LDB );
    }

    free( A );
    free( B );

    return hres;
}

testing_t   test_zlatro;
const char *zlatro_params[] = { "mtxfmt", "nb", "trans", "uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zlatro_output[] = { NULL };
const char *zlatro_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlatro_init( void ) __attribute__( ( constructor ) );
void
testing_zlatro_init( void )
{
    test_zlatro.name   = "zlatro";
    test_zlatro.helper = "Matrix Transposition";
    test_zlatro.params = zlatro_params;
    test_zlatro.output = zlatro_output;
    test_zlatro.outchk = zlatro_outchk;
    test_zlatro.fptr_desc = testing_zlatro_desc;
    test_zlatro.fptr_std  = testing_zlatro_std;
    test_zlatro.next   = NULL;

    testing_register( &test_zlatro );
}
