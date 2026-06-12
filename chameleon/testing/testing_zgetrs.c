/**
 *
 * @file testing_zgetrs.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf testing
 *
 * @version 1.4.0
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "chameleon/chameleon_z.h"
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <chameleon/getenv.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/lapacke.h>
#endif

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zgetrs_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;
    int        P, Q;

    /* Read arguments */
    int           async = parameters_getvalue_int( "async" );
    int           nb    = run_arg_get_nb( args );
    int           ib    = run_arg_get_ib( args );
    cham_trans_t  trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    int           N     = run_arg_get_int( args, "N", 1000 );
    int           NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int           LDA   = run_arg_get_int( args, "LDA", N );
    int           LDB   = run_arg_get_int( args, "LDB", N );
    int           seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int           seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;
    CHAM_ipiv_t *descIPIV;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );
    if ( descA->get_rankof_init != chameleon_getrankof_2d ) {
        parameters_desc_destroy( &descA );
        fprintf( stderr,
                 "SKIPPED: zgetrf/zlaswp kernels only supports 2DBC data distributions\n" );
        return -1;
    }

    parameters_desc_create( "X", &descX, ChamComplexDouble, nb, nb, LDB, NRHS, N, NRHS );

    P = chameleon_desc_datadist_get_iparam( descA, 0 );
    Q = chameleon_desc_datadist_get_iparam( descA, 1 );

    CHAMELEON_Ipiv_Create( &descIPIV, ChamLeft, descA->mb, N, P, P*Q, NULL );

    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    CHAMELEON_zgetrf_Tile( descA, descIPIV );

    if ( async ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( ChamLeft, descA );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgetrs_Tile_Async( trans, descA, descIPIV, descX, ws, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgetrs_Tile( trans, descA, descIPIV, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ws != NULL ) {
        CHAMELEON_zlaswp_WS_Free( ws );
    }

    /* Checks the factorization and residual */
#if !defined(CHAMELEON_SIMULATION)
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0, *descB;

        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        descB  = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB,  seedB );

        hres += check_zsolve( args, ChamGeneral, trans, ChamUpperLower, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    CHAMELEON_Ipiv_Destroy( &descIPIV );
    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}
#endif

int
testing_zgetrs_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int           api   = parameters_getvalue_int( "api" );
#endif
    int           nb    = run_arg_get_nb( args );
    cham_trans_t  trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    int           N     = run_arg_get_int( args, "N", 1000 );
    int           NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int           LDA   = run_arg_get_int( args, "LDA", N );
    int           LDB   = run_arg_get_int( args, "LDB", N );
    int           seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int           seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;
    int *IPIV;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );
    IPIV = malloc( sizeof(int) * N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( N, N, A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    CHAMELEON_zgetrf( N, N, A, LDA, IPIV );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zgetrs( LAPACK_COL_MAJOR, chameleon_lapack_const(trans), N, NRHS, A, LDA, IPIV, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zgetrs( trans, N, NRHS, A, LDA, IPIV, X, LDB);
        break;
#if !defined(CHAMELEON_SIMULATION) & 0
    case 2:
        CHAMELEON_lapacke_zgetrs( CblasColMajor, N, NRHS, A, LDA, IPIV, B, LDB );
        break;
#endif
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrs( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

#if !defined(CHAMELEON_SIMULATION)
    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        CHAMELEON_zplrnt( N, N, A0, LDA, seedA );
        CHAMELEON_zplrnt( N, NRHS, B,  LDB, seedB );

        hres += check_zsolve_std( args, ChamGeneral, trans, ChamUpperLower, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }
#endif
#endif

    free ( IPIV );
    free( A );
    free( X );

    (void)check;
    return hres;
}

testing_t   test_zgetrs;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zgetrs_params[] = { "m", "n", "lda", "seedA", NULL };
#else
const char *zgetrs_params[] = { "mtxfmt", "nb", "ib", "trans", "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
#endif
const char *zgetrs_output[] = { NULL };
const char *zgetrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgetrs_init( void ) __attribute__( ( constructor ) );
void
testing_zgetrs_init( void )
{
    test_zgetrs.name   = "zgetrs";
    test_zgetrs.helper = "General triangular solve (LU with partial pivoting)";
    test_zgetrs.params = zgetrs_params;
    test_zgetrs.output = zgetrs_output;
    test_zgetrs.outchk = zgetrs_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zgetrs.fptr_desc = NULL;
#else
    test_zgetrs.fptr_desc = testing_zgetrs_desc;
#endif
    test_zgetrs.fptr_std  = testing_zgetrs_std;
    test_zgetrs.next   = NULL;

    testing_register( &test_zgetrs );
}
