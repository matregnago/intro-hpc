/**
 *
 * @file testing_zgeqrs.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrs testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
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

int
testing_zgeqrs_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      LDB    = run_arg_get_int( args, "LDB", M );
    int      RH     = run_arg_get_int( args, "qra", 0 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      seedB  = run_arg_get_int( args, "seedB", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descX, *descT;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( N > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: The QR solution is performed only when M >= N\n" );
        }
        return -1;
    }

    if ( RH > 0 ) {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamTreeHouseholder );
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_SIZE, RH );
    }
    else {
        CHAMELEON_Set( CHAMELEON_HOUSEHOLDER_MODE, ChamFlatHouseholder );
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descX, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDB, NRHS, 0, 0, M, NRHS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &descT, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    hres = CHAMELEON_zgeqrf_Tile( descA, descT );

    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgeqrs_Tile_Async( descA, descT, descX,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgeqrs_Tile( descA, descT, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgeqrs( M, N, NRHS ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation, orthogonality and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *descB  = CHAMELEON_Desc_Copy( descX, CHAMELEON_MAT_ALLOC_TILE );
        CHAM_desc_t *subX   = CHAMELEON_Desc_SubMatrix( descX, 0, 0, N, NRHS );
        CHAM_desc_t *subB   = CHAMELEON_Desc_SubMatrix( descB, 0, 0, M, NRHS );

        CHAMELEON_zplrnt_Tile( descA0, seedA );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        hres += check_zsolve( args, ChamGeneral, ChamNoTrans, ChamUpperLower, descA0, subX, subB );

        free( subB );
        free( subX );
        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

testing_t   test_zgeqrs;
const char *zgeqrs_params[] = { "mtxfmt", "nb",  "ib",  "m",     "n",     "k",
                                "lda",    "ldb", "qra", "seedA", "seedB", NULL };
const char *zgeqrs_output[] = { NULL };
const char *zgeqrs_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgeqrs_init( void ) __attribute__( ( constructor ) );
void
testing_zgeqrs_init( void )
{
    test_zgeqrs.name   = "zgeqrs";
    test_zgeqrs.helper = "General QR solve";
    test_zgeqrs.params = zgeqrs_params;
    test_zgeqrs.output = zgeqrs_output;
    test_zgeqrs.outchk = zgeqrs_outchk;
    test_zgeqrs.fptr_desc = testing_zgeqrs_desc;
    test_zgeqrs.fptr_std  = NULL;
    test_zgeqrs.next   = NULL;

    testing_register( &test_zgeqrs );
}
