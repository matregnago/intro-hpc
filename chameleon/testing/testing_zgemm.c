/**
 *
 * @file testing_zgemm.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm testing
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
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zgemm_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    int          nb     = run_arg_get_nb( args );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t transB = run_arg_get_trans( args, "transB", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( transA == ChamNoTrans ) ? M : K ) );
    int          LDB    = run_arg_get_int( args, "LDB", ( ( transB == ChamNoTrans ) ? K : N ) );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An, Bm, Bn;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;
    void        *ws = NULL;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Am = M;
        An = K;
    }
    else {
        Am = K;
        An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K;
        Bn = N;
    }
    else {
        Bm = N;
        Bn = K;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, An, Am, An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, Bn, Bm, Bn );
    parameters_desc_create( "C", &descC, ChamComplexDouble, nb, nb, LDC, N, M, N );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    if ( async ) {
        ws = CHAMELEON_zgemm_WS_Alloc( transA, transB, descA, descB, descC );
    }

    /* Calculate the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgemm_Tile_Async( transA, transB, alpha, descA, descB, beta, descC, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgemm_Tile( transA, transB, alpha, descA, descB, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descCinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDC, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres += check_zgemm( args, transA, transB, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );
    parameters_desc_destroy( &descC );

    return hres;
}
#endif

int
testing_zgemm_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int          api    = parameters_getvalue_int( "api" );
#endif
    int          nb     = run_arg_get_nb( args );
    cham_trans_t transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t transB = run_arg_get_trans( args, "transB", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( transA == ChamNoTrans ) ? M : K ) );
    int          LDB    = run_arg_get_int( args, "LDB", ( ( transB == ChamNoTrans ) ? K : N ) );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   seedC = run_arg_get_int( args, "seedC", testing_ialea() );

    /* Descriptors */
    int                    Am, An, Bm, Bn;
    CHAMELEON_Complex64_t *A, *B, *C;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Am = M;
        An = K;
    }
    else {
        Am = K;
        An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K;
        Bn = N;
    }
    else {
        Bm = N;
        Bn = K;
    }

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*Bn );
    C = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplrnt( Bm, Bn, B, LDB, seedB );
    CHAMELEON_zplrnt( M, N, C, LDC, seedC );

    /* Calculate the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                        CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), C, LDC );
    testing_stop( &test_data, flops_zgemm( M, N, K ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zgemm( transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
        break;
    case 2:
        CHAMELEON_cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                               CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), C, LDC );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Cinit;
        Cinit = malloc( sizeof(CHAMELEON_Complex64_t) * LDC*N );
        CHAMELEON_zplrnt( M, N, Cinit, LDC, seedC );

        hres += check_zgemm_std( args, transA, transB, alpha, M, N, K, A, LDA, B, LDB, beta, Cinit, C, LDC );

        free( Cinit );
    }
#endif

    free( A );
    free( B );
    free( C );

    (void)check;
    return hres;
}

testing_t   test_zgemm;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zgemm_params[] = { "transA", "transB", "m",     "n",
                               "k",      "lda",   "ldb",    "ldc",    "alpha", "beta",
                               "seedA",  "seedB", "seedC",  NULL };
#else
const char *zgemm_params[] = { "mtxfmt", "nb",    "transA", "transB", "m",     "n",
                               "k",      "lda",   "ldb",    "ldc",    "alpha", "beta",
                               "seedA",  "seedB", "seedC",  NULL };
#endif
const char *zgemm_output[] = { NULL };
const char *zgemm_outchk[] = { "||A||", "||B||", "||C||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgemm_init( void ) __attribute__( ( constructor ) );
void
testing_zgemm_init( void )
{
    test_zgemm.name   = "zgemm";
    test_zgemm.helper = "General matrix-matrix multiply";
    test_zgemm.params = zgemm_params;
    test_zgemm.output = zgemm_output;
    test_zgemm.outchk = zgemm_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zgemm.fptr_desc = NULL;
#else
    test_zgemm.fptr_desc = testing_zgemm_desc;
#endif
    test_zgemm.fptr_std  = testing_zgemm_std;
    test_zgemm.next   = NULL;

    testing_register( &test_zgemm );
}
