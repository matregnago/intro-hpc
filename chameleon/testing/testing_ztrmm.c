/**
 *
 * @file testing_ztrmm.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm testing
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
#include <coreblas/cblas.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_ztrmm_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int                   async = parameters_getvalue_int( "async" );
    int                   nb    = run_arg_get_nb( args );
    int                   P     = parameters_getvalue_int( "P" );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( side == ChamLeft ) ? M : N );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descB, *descBinit;

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the side */
    if ( side == ChamLeft ) {
        parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, M, M, M );
    }
    else {
        parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );
    }

    /* Creates the matrices */
    parameters_desc_create( "B", &descB, ChamComplexDouble, nb, nb, LDB, N, M, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Calculates the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_ztrmm_Tile_Async( side, uplo, trans, diag, alpha, descA, descB,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_ztrmm_Tile( side, uplo, trans, diag, alpha, descA, descB );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrmm( side, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Desc_Create(
            &descBinit, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, LDB, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descBinit, seedB );

        hres += check_ztrmm( args, CHECK_TRMM, side, uplo, trans, diag,
                             alpha, descA, descB, descBinit );

        CHAMELEON_Desc_Destroy( &descBinit );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descB );

    return hres;
}
#endif

int
testing_ztrmm_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int                   api   = parameters_getvalue_int( "api" );
#endif
    int                   nb    = run_arg_get_nb( args );
    cham_trans_t          trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t           side  = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t           uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t           diag  = run_arg_get_diag( args, "diag", ChamNonUnit );
    int                   N     = run_arg_get_int( args, "N", 1000 );
    int                   M     = run_arg_get_int( args, "M", N );
    int                   LDA   = run_arg_get_int( args, "LDA", ( side == ChamLeft ) ? M : N );
    int                   LDB   = run_arg_get_int( args, "LDB", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int                   seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    int                    An;
    CHAMELEON_Complex64_t *A, *B;

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculates the dimensions according to the side */
    An = ( side == ChamLeft ) ? M : N;

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*An );
    B = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N  );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( An, An, A, LDA, seedA );
    CHAMELEON_zplrnt( M,  N,  B, LDB, seedB );

    /* Calculates the product */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    cblas_ztrmm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 (CBLAS_DIAG)diag, M, N, CBLAS_SADDR(alpha), A, LDA, B, LDB );
    testing_stop( &test_data, flops_ztrmm( side, M, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_ztrmm( side, uplo, trans, diag, M, N, alpha, A, LDA, B, LDB );
        break;
    case 2:
        CHAMELEON_cblas_ztrmm( CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                               (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag, M, N,
                               CBLAS_SADDR(alpha), A, LDA, B, LDB );
        break;
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_ztrmm( side, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *Binit;
        Binit = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*N );
        CHAMELEON_zplrnt( M, N, Binit, LDB, seedB );

        hres += check_ztrmm_std( args, CHECK_TRMM, side, uplo, trans, diag, M, N, alpha, A, LDA, B, Binit, LDB );

        free( Binit );
    }
#endif

    free( A );
    free( B );

    (void)check;
    return hres;
}

testing_t   test_ztrmm;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *ztrmm_params[] = { "trans", "side",  "uplo",  "diag",  "m",
                               "n",      "lda", "ldb",   "alpha", "seedA", "seedB", NULL };
#else
const char *ztrmm_params[] = { "mtxfmt", "nb",  "trans", "side",  "uplo",  "diag",  "m",
                               "n",      "lda", "ldb",   "alpha", "seedA", "seedB", NULL };
#endif
const char *ztrmm_output[] = { NULL };
const char *ztrmm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrmm_init( void ) __attribute__( ( constructor ) );
void
testing_ztrmm_init( void )
{
    test_ztrmm.name   = "ztrmm";
    test_ztrmm.helper = "Triangular matrix-matrix multiply";
    test_ztrmm.params = ztrmm_params;
    test_ztrmm.output = ztrmm_output;
    test_ztrmm.outchk = ztrmm_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_ztrmm.fptr_desc = NULL;
#else
    test_ztrmm.fptr_desc = testing_ztrmm_desc;
#endif
    test_ztrmm.fptr_std  = testing_ztrmm_std;
    test_ztrmm.next   = NULL;

    testing_register( &test_ztrmm );
}
