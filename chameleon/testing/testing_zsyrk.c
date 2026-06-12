/**
 *
 * @file testing_zsyrk.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> z c d s
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
testing_zsyrk_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    int                   P     = parameters_getvalue_int( "P" );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   K     = run_arg_get_int( args, "K", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    CHAMELEON_Complex64_t bump  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC, *descCinit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );
    bump  = run_arg_get_complex64( args, "bump", bump );

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
    CHAMELEON_zplgsy_Tile( bump, uplo, descC, seedC );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zsyrk_Tile_Async( uplo, trans, alpha, descA, beta, descC,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zsyrk_Tile( uplo, trans, alpha, descA, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsyrk( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descCinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplgsy_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamSymmetric, uplo, trans, alpha, descA, NULL,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descC );

    return hres;
}
#endif

int
testing_zsyrk_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int          api   = parameters_getvalue_int( "api" );
#endif
    int          nb    = run_arg_get_nb( args );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          K     = run_arg_get_int( args, "K", N );
    int          LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int          LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    CHAMELEON_Complex64_t bump  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );

    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *C;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );
    bump  = run_arg_get_complex64( args, "bump", bump );

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
    CHAMELEON_zplgsy( bump, uplo, N, C, LDC, seedC );

    /* Calculates the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_zsyrk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                 CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), C, LDC );
    testing_stop( &test_data, flops_zsyrk( K, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zsyrk( uplo, trans, N, K, alpha, A, LDA, beta, C, LDC );
        break;
    case 2:
        CHAMELEON_cblas_zsyrk( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                               CBLAS_SADDR(alpha), A, LDA, CBLAS_SADDR(beta), C, LDC );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zsyrk( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Cinit;
        Cinit = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );
        CHAMELEON_zplgsy( bump, uplo, N, Cinit, LDC, seedC );

        hres += check_zsyrk_std( args, ChamSymmetric, uplo, trans, N, K, alpha, A, LDA, NULL, LDA, beta, Cinit, C, LDC );

        free( Cinit );
    }
#endif

    free( A );
    free( C );

    (void)check;
    return hres;
}

testing_t   test_zsyrk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zsyrk_params[] = { "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
#else
const char *zsyrk_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",     "k",    "lda",
                               "ldc",    "alpha", "beta",  "seedA", "seedC", "bump", NULL };
#endif
const char *zsyrk_output[] = { NULL };
const char *zsyrk_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsyrk_init( void ) __attribute__( ( constructor ) );
void
testing_zsyrk_init( void )
{
    test_zsyrk.name   = "zsyrk";
    test_zsyrk.helper = "Symmetrix matrix-matrix rank k update";
    test_zsyrk.params = zsyrk_params;
    test_zsyrk.output = zsyrk_output;
    test_zsyrk.outchk = zsyrk_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zsyrk.fptr_desc = NULL;
#else
    test_zsyrk.fptr_desc = testing_zsyrk_desc;
#endif
    test_zsyrk.fptr_std  = testing_zsyrk_std;
    test_zsyrk.next   = NULL;

    testing_register( &test_zsyrk );
}
