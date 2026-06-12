/**
 *
 * @file testing_zgeadd.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgeadd( int _M, int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t M     = _M;
    cham_fixdbl_t N     = _N;

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* 2 multiplications and 1 addition per element */
    flops = ( 2. * 6. + 2. ) * M * N;
#else
    flops = ( 2. + 1. ) * M * N;
#endif

    return flops;
}

int
testing_zgeadd_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    if ( trans != ChamNoTrans ) {
        Am = N;
        An = M;
    }
    else {
        Am = M;
        An = N;
    }

    /* Create the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, An, Am, An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, N, M, N );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Compute the sum */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgeadd_Tile_Async( trans, alpha, descA, beta, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgeadd_Tile( trans, alpha, descA, beta, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeadd( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descB0 = CHAMELEON_Desc_Copy( descB, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descB0, seedB );

        hres += check_zsum( args, ChamUpperLower, trans, alpha, descA, beta, descB0, descB );

        CHAMELEON_Desc_Destroy( &descB0 );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}

int
testing_zgeadd_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int          nb     = run_arg_get_nb( args );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int          LDB    = run_arg_get_int( args, "LDB", M );
    int          seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedB  = run_arg_get_int( args, "seedB", testing_ialea() );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();

    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *B;

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    Am = (trans == ChamNoTrans)? M : N;
    An = (trans == ChamNoTrans)? N : M;

    /* Create the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N  );

    /* Fill the matrix with random values */
    CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
    CHAMELEON_zplrnt( M,  N,  B, LDB, seedB );

    /* Compute the sum */
    testing_start( &test_data );
    hres = CHAMELEON_zgeadd( trans, M, N, alpha, A, LDA, beta, B, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeadd( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Check the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *B0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N  );
        CHAMELEON_zplrnt( M, N, B0, LDB, seedB );

        hres += check_zsum_std( args, ChamUpperLower, trans, M, N, alpha, A, LDA, beta, B0, B, LDB );

        free( B0 );
    }

    free( A );
    free( B );

    return hres;
}

testing_t   test_zgeadd;
const char *zgeadd_params[] = { "mtxfmt", "nb",    "trans", "m",     "n",     "lda",
                                "ldb",    "alpha", "beta",  "seedA", "seedB", NULL };
const char *zgeadd_output[] = { NULL };
const char *zgeadd_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeadd_init( void ) __attribute__( ( constructor ) );
void
testing_zgeadd_init( void )
{
    test_zgeadd.name   = "zgeadd";
    test_zgeadd.helper = "General matrix-matrix addition";
    test_zgeadd.params = zgeadd_params;
    test_zgeadd.output = zgeadd_output;
    test_zgeadd.outchk = zgeadd_outchk;
    test_zgeadd.fptr_desc = testing_zgeadd_desc;
    test_zgeadd.fptr_std  = testing_zgeadd_std;
    test_zgeadd.next   = NULL;

    testing_register( &test_zgeadd );
}
