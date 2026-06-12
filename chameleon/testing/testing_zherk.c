/**
 *
 * @file testing_zherk.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> z c
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/cblas.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zherk_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async  = parameters_getvalue_int( "async" );
    int                   nb     = run_arg_get_nb( args );
    int                   P      = parameters_getvalue_int( "P" );
    cham_trans_t          trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N      = run_arg_get_int( args, "N", 1000 );
    int                   K      = run_arg_get_int( args, "K", N );
    int                   LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC    = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t zalpha = testing_zalea();
    CHAMELEON_Complex64_t zbeta  = testing_zalea();
    CHAMELEON_Complex64_t zbump  = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedC  = run_arg_get_int( args, "seedC", testing_ialea() );
    int                   Q      = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC, *descCinit;
    double       alpha, beta, bump;

    zalpha = run_arg_get_complex64( args, "alpha", zalpha );
    zbeta  = run_arg_get_complex64( args, "beta",  zbeta );
    zbump  = run_arg_get_complex64( args, "bump",  zbump );
    alpha  = creal( zalpha );
    beta   = creal( zbeta );
    bump   = creal( zbump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, An, Am, An );
    parameters_desc_create( "C", &descC, ChamComplexDouble, nb, nb, LDC, N, N, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplghe_Tile( bump, uplo, descC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zherk_Tile_Async( uplo, trans, alpha, descA, beta, descC,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zherk_Tile( uplo, trans, alpha, descA, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zherk( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descCinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, NULL,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descC );

    return hres;
}
#endif

int
testing_zherk_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int                   api    = parameters_getvalue_int( "api" );
#endif
    int                   nb     = run_arg_get_nb( args );
    cham_trans_t          trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N      = run_arg_get_int( args, "N", 1000 );
    int                   K      = run_arg_get_int( args, "K", N );
    int                   LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC    = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t zalpha = testing_zalea();
    CHAMELEON_Complex64_t zbeta  = testing_zalea();
    CHAMELEON_Complex64_t zbump  = testing_zalea();
    int                   seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedC  = run_arg_get_int( args, "seedC", testing_ialea() );


    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *C;
    double                 alpha, beta, bump;

    zalpha = run_arg_get_complex64( args, "alpha", zalpha );
    zbeta  = run_arg_get_complex64( args, "beta",  zbeta );
    zbump  = run_arg_get_complex64( args, "bump",  zbump );
    alpha  = creal( zalpha );
    beta   = creal( zbeta );
    bump   = creal( zbump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    C = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N  );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplghe( bump, uplo, N, C, LDC, seedC );

    /* Calculates the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_zherk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                 alpha, A, LDA, beta, C, LDC );
    testing_stop( &test_data, flops_zherk( K, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zherk( uplo, trans, N, K, alpha, A, LDA, beta, C, LDC );
        break;
    case 2:
        CHAMELEON_cblas_zherk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                               alpha, A, LDA, beta, C, LDC );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zherk( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Cinit;
        Cinit = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );
        CHAMELEON_zplghe( bump, uplo, N, Cinit, LDC, seedC );

        hres += check_zsyrk_std( args, ChamHermitian, uplo, trans, N, K, alpha, A, LDA, NULL, LDA, beta, Cinit, C, LDC );

        free( Cinit );
    }
#endif

    free( A );
    free( C );

    (void)check;
    return hres;
}

testing_t   test_zherk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zherk_params[] = { "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
#else
const char *zherk_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
#endif
const char *zherk_output[] = { NULL };
const char *zherk_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zherk_init( void ) __attribute__( ( constructor ) );
void
testing_zherk_init( void )
{
    test_zherk.name   = "zherk";
    test_zherk.helper = "Hermitian matrix-matrix rank k update";
    test_zherk.params = zherk_params;
    test_zherk.output = zherk_output;
    test_zherk.outchk = zherk_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zherk.fptr_desc = NULL;
#else
    test_zherk.fptr_desc = testing_zherk_desc;
#endif
    test_zherk.fptr_std  = testing_zherk_std;
    test_zherk.next   = NULL;

    testing_register( &test_zherk );
}
