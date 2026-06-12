/**
 *
 * @file testing_zher2k.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k testing
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
testing_zher2k_desc( run_arg_list_t *args, int check )
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
    int                   LDB   = run_arg_get_int( args, "LDB", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t zbeta = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    CHAMELEON_Complex64_t zbump = testing_zalea();
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;
    double       beta, bump;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    zbeta = run_arg_get_complex64( args, "beta", zbeta );
    zbump = run_arg_get_complex64( args, "bump", zbump );
    beta  = creal( zbeta );
    bump  = creal( zbump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, An, Am, An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, An, Am, An );
    parameters_desc_create( "C", &descC, ChamComplexDouble, nb, nb, LDC, N, N, N );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplghe_Tile( bump, uplo, descC, seedC );

    /* Calculate the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zher2k_Tile_Async( uplo, trans, alpha, descA, descB, beta, descC,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zher2k_Tile( uplo, trans, alpha, descA, descB, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zher2k( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descCinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( bump, uplo, descCinit, seedC );

        hres += check_zsyrk( args, ChamHermitian, uplo, trans, alpha, descA, descB,
                             beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );
    parameters_desc_destroy( &descC );

    return hres;
}
#endif

int
testing_zher2k_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int                   api   = parameters_getvalue_int( "api" );
#endif
    int                   nb    = run_arg_get_nb( args );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   K     = run_arg_get_int( args, "K", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDB   = run_arg_get_int( args, "LDB", ( ( trans == ChamNoTrans ) ? N : K ) );
    int                   LDC   = run_arg_get_int( args, "LDC", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t zbeta = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    CHAMELEON_Complex64_t zbump = testing_zalea();

    /* Descriptors */
    int                    Am, An, Bm, Bn;
    CHAMELEON_Complex64_t *A, *B, *C;
    double                 beta, bump;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    zbeta = run_arg_get_complex64( args, "beta", zbeta );
    zbump = run_arg_get_complex64( args, "bump", zbump );
    beta  = creal( zbeta );
    bump  = creal( zbump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
        Bm = N;
        Bn = K;
    }
    else {
        Am = K;
        An = N;
        Bm = K;
        Bn = N;
    }

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*Bn );
    C = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N  );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplrnt( Bm, Bn, B, LDB, seedB );
    CHAMELEON_zplghe( bump, uplo, N, C, LDC, seedC );

    /* Calculate the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_zher2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                  CBLAS_SADDR(alpha), A, LDA, B, LDB, beta, C, LDC );
    testing_stop( &test_data, flops_zher2k( K, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zher2k( uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
        break;
    case 2:
        CHAMELEON_cblas_zher2k( CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, N, K,
                                CBLAS_SADDR(alpha), A, LDA, B, LDB, beta, C, LDC );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zher2k( K, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Cinit;
        Cinit = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );
        CHAMELEON_zplghe( bump, uplo, N, Cinit, LDC, seedC );

        hres += check_zsyrk_std( args, ChamHermitian, uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, Cinit, C, LDC );

        free( Cinit );
    }
#endif

    free( A );
    free( B );
    free( C );

    (void)check;
    return hres;
}

testing_t   test_zher2k;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zher2k_params[] = { "trans", "uplo",  "n",    "k",
                                "lda",    "ldb",   "ldc",   "alpha", "beta", "seedA",
                                "seedB",  "seedC", "bump",  NULL };
#else
const char *zher2k_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "n",    "k",
                                "lda",    "ldb",   "ldc",   "alpha", "beta", "seedA",
                                "seedB",  "seedC", "bump",  NULL };
#endif
const char *zher2k_output[] = { NULL };
const char *zher2k_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zher2k_init( void ) __attribute__( ( constructor ) );
void
testing_zher2k_init( void )
{
    test_zher2k.name   = "zher2k";
    test_zher2k.helper = "Hermitian matrix-matrix rank 2k update";
    test_zher2k.params = zher2k_params;
    test_zher2k.output = zher2k_output;
    test_zher2k.outchk = zher2k_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zher2k.fptr_desc = NULL;
#else
    test_zher2k.fptr_desc = testing_zher2k_desc;
#endif
    test_zher2k.fptr_std  = testing_zher2k_std;
    test_zher2k.next   = NULL;

    testing_register( &test_zher2k );
}
