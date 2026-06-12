/**
 *
 * @file testing_zgesv.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesv testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "testings.h"
#include "chameleon/chameleon_z.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <chameleon/getenv.h>
#include <coreblas/lapacke.h>

static cham_fixdbl_t
flops_zgesv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zgetrf( N, N ) + flops_zgetrs( N, NRHS );
    return flops;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zgesv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;
    int        P, Q;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    int         ib    = run_arg_get_ib( args );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         LDB   = run_arg_get_int( args, "LDB", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int         seedB = run_arg_get_int( args, "seedB", testing_ialea() );

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

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    if ( async ) {
        ws = CHAMELEON_zgetrf_WS_Alloc( descA );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgesv_Tile_Async( descA, descIPIV, descX, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgesv_Tile( descA, descIPIV, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgesv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( async ) {
        CHAMELEON_zgetrf_WS_Free( ws );
    }

    /* Checks the factorisation and the residual */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        CHAMELEON_zlaswp_Tile( ChamLeft, ChamDirForward, descA0, 1, N, descIPIV );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        if ( hres ) {
            CHAMELEON_Desc_Destroy( &descA0 );
            CHAMELEON_Ipiv_Destroy( &descIPIV );
            parameters_desc_destroy( &descA );
            parameters_desc_destroy( &descX );
            return hres;
        }

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        hres += check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Ipiv_Destroy( &descIPIV );
    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descX );

    return hres;
}
#endif

int
testing_zgesv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int    api   = parameters_getvalue_int( "api" );
#endif
    int    nb    = run_arg_get_nb( args );
    int    ib    = run_arg_get_ib( args );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    NRHS  = run_arg_get_int( args, "NRHS", 1 );
    int    LDA   = run_arg_get_int( args, "LDA", N );
    int    LDB   = run_arg_get_int( args, "LDB", N );
    int    seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int    seedB = run_arg_get_int( args, "seedB", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *X;
    int *IPIV;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    X = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );
    IPIV = malloc( sizeof(int) * N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( N, N,    A, LDA, seedA );
    CHAMELEON_zplrnt( N, NRHS, X, LDB, seedB );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zgesv( LAPACK_COL_MAJOR, N, NRHS, A, LDA, IPIV, X, LDB );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgesv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zgesv( N, NRHS, A, LDA, IPIV, X, LDB );
        break;
#if !defined(CHAMELEON_SIMULATION) && 0
    case 2:
        CHAMELEON_lapacke_zgesv( CblasColMajor, chameleon_lapack_const(uplo), N, NRHS, A, LDA, X, LDB );
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
    testing_stop( &test_data, flops_zgesv( N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residual */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAMELEON_Complex64_t *A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_Complex64_t *B  = malloc( sizeof(CHAMELEON_Complex64_t) * LDB*NRHS );

        /* Check the factorization */
        CHAMELEON_zplrnt( N, N, A0, LDA, seedA );
        CHAMELEON_zlaswp( ChamLeft, ChamDirForward, N, N, A0, LDA, 1, N, IPIV );

        hres += check_zxxtrf_std( args, ChamGeneral, ChamUpperLower, N, N, A0, A, LDA );

        /* Check the solve */
        CHAMELEON_zplrnt( N, N,    A0, LDA, seedA );
        CHAMELEON_zplrnt( N, NRHS, B,  LDB, seedB );
        hres += check_zsolve_std( args, ChamGeneral, ChamNoTrans, ChamUpperLower, N, NRHS, A0, LDA, X, B, LDB );

        free( A0 );
        free( B );
    }
#endif

    free( A );
    free( X );

    (void)check;
    return hres;
}

testing_t   test_zgesv;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zgesv_params[] = { "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
#else
const char *zgesv_params[] = { "mtxfmt", "nb", "ib", "n", "nrhs", "lda", "ldb", "seedA", "seedB", NULL };
#endif
const char *zgesv_output[] = { NULL };
const char *zgesv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgesv_init( void ) __attribute__( ( constructor ) );
void
testing_zgesv_init( void )
{
    test_zgesv.name      = "zgesv";
    test_zgesv.helper    = "General linear system solve (LU with partial pivoting)";
    test_zgesv.params    = zgesv_params;
    test_zgesv.output    = zgesv_output;
    test_zgesv.outchk    = zgesv_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zgesv.fptr_desc = NULL;
#else
    test_zgesv.fptr_desc = testing_zgesv_desc;
#endif
    test_zgesv.fptr_std  = testing_zgesv_std;
    test_zgesv.next      = NULL;

    testing_register( &test_zgesv );
}

