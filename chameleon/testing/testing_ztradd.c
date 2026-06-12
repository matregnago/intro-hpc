/**
 *
 * @file testing_ztradd.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd testing
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
flops_ztradd( cham_uplo_t uplo, int _M, int _N )
{
    cham_fixdbl_t flops = 0.;
    cham_fixdbl_t minMN = (cham_fixdbl_t)chameleon_min( _M, _N );
    cham_fixdbl_t M     = _M;
    cham_fixdbl_t N     = _N;

    switch ( uplo ) {
        case ChamUpper:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + M * (cham_fixdbl_t)chameleon_max( 0, _N - _M );
            break;
        case ChamLower:
            flops = ( minMN * ( minMN + 1 ) / 2 ) + N * (cham_fixdbl_t)chameleon_max( 0, _M - _N );
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    /* 2 multiplications and 1 addition per element */
    flops *= ( 2. * 6. + 2. );
#else
    flops *= ( 2. + 1. );
#endif

    return flops;
}

int
testing_ztradd_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descB;
    cham_uplo_t  uplo_inv = uplo;

    if ( uplo != ChamUpperLower && trans != ChamNoTrans ) {
        uplo_inv = (uplo == ChamUpper) ? ChamLower : ChamUpper;
    }

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

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, An, Am, An );
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, N, M, N );

    /* Fills the matrix with random values */
    switch ( uplo ) {
        case ChamUpper:
        case ChamLower:
            CHAMELEON_zplgsy_Tile( 0., uplo_inv, descA, seedA );
            CHAMELEON_zplgsy_Tile( 0., uplo,     descB, seedB );
            break;
        case ChamUpperLower:
        default:
            CHAMELEON_zplrnt_Tile( descA, seedA );
            CHAMELEON_zplrnt_Tile( descB, seedB );
            break;
    }

    /* Calculates the sum */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_ztradd_Tile_Async( uplo, trans, alpha, descA, beta, descB,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_ztradd_Tile( uplo, trans, alpha, descA, beta, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztradd( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descB0 = CHAMELEON_Desc_Copy( descB, CHAMELEON_MAT_ALLOC_TILE );

        if ( uplo == ChamUpperLower ) {
            CHAMELEON_zplrnt_Tile( descB0, seedB );
        }
        else {
            CHAMELEON_zplgsy_Tile( 0., uplo, descB0, seedB );
        }
        hres += check_zsum( args, uplo, trans, alpha, descA, beta, descB0, descB );

        CHAMELEON_Desc_Destroy( &descB0 );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}

int
testing_ztradd_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int          nb    = run_arg_get_nb( args );
    cham_trans_t trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t  uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int          N     = run_arg_get_int( args, "N", 1000 );
    int          M     = run_arg_get_int( args, "M", N );
    int          LDA   = run_arg_get_int( args, "LDA", ( ( trans == ChamNoTrans ) ? M : N ) );
    int          LDB   = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int          seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int          seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    int                    Am, An;
    CHAMELEON_Complex64_t *A, *B;
    cham_uplo_t            uplo_inv = uplo;

    if ( (uplo != ChamUpperLower) && (trans != ChamNoTrans) ) {
        uplo_inv = (uplo == ChamUpper) ? ChamLower : ChamUpper;
    }

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

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N );

    /* Fills the matrix with random values */
    switch ( uplo ) {
        case ChamUpper:
        case ChamLower:
            CHAMELEON_zplgtr( 0., uplo_inv, Am, An, A, LDA, seedA );
            CHAMELEON_zplgtr( 0., uplo,     M,  N,  B, LDB, seedB );
            break;
        case ChamUpperLower:
        default:
            CHAMELEON_zplrnt( Am, An, A, LDA, seedA );
            CHAMELEON_zplrnt( M,  N,  B, LDB, seedB );
            break;
    }

    /* Calculates the sum */
    testing_start( &test_data );
    hres = CHAMELEON_ztradd( uplo, trans, M, N, alpha, A, LDA, beta, B, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztradd( uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *B0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N );

        if ( uplo == ChamUpperLower ) {
            CHAMELEON_zplrnt( M, N, B0, LDB, seedB );
        }
        else {
            CHAMELEON_zplgtr( 0., uplo, M, N, B0, LDB, seedB );
        }
        hres += check_zsum_std( args, uplo, trans, M, N, alpha, A, LDA, beta, B0, B, LDB );

        free( B0 );
    }

    free( A );
    free( B );

    return hres;
}

testing_t   test_ztradd;
const char *ztradd_params[] = { "mtxfmt", "nb",    "trans", "uplo",  "m",     "n", "lda",
                                "ldb",    "alpha", "beta",  "seedA", "seedB", NULL };
const char *ztradd_output[] = { NULL };
const char *ztradd_outchk[] = {  "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztradd_init( void ) __attribute__( ( constructor ) );
void
testing_ztradd_init( void )
{
    test_ztradd.name   = "ztradd";
    test_ztradd.helper = "Triangular matrix-matrix addition";
    test_ztradd.params = ztradd_params;
    test_ztradd.output = ztradd_output;
    test_ztradd.outchk = ztradd_outchk;
    test_ztradd.fptr_desc = testing_ztradd_desc;
    test_ztradd.fptr_std  = testing_ztradd_std;
    test_ztradd.next   = NULL;

    testing_register( &test_ztradd );
}
