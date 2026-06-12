/**
 *
 * @file testing_zgepdf_qr.c
 *
 * @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2020-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qr testing
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @author Alycia Lisito
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include <chameleon/getenv.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#endif
#include <libhqr.h>

static cham_fixdbl_t
flops_zgepdf_qr( int M, int N )
{
    double flops = flops_zgeqrf( M + N, N ) + flops_zungqr( M + N, N, N );
    return flops;
}

int
testing_zgepdf_qr_desc( run_arg_list_t *args, int check )
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
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t  *descA1, *descA2, *descQ1, *descQ2;
    CHAM_desc_t  *TS1, *TT1, *TS2, *TT2;
    libhqr_tree_t qrtreeT, qrtreeB;
    int           zqdwh_opt_id = chameleon_env_on_off( "CHAMELEON_TESTING_GEPDF_OPTID", CHAMELEON_TRUE );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( N > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The QR factorization for Polar Decomposition is performed only when "
                     "M >= N\n" );
        }
        return -1;
    }

    if ( ( N % nb ) != 0 ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The QR factorization for Polar Decomposition supports only multiple "
                     "of nb\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA1, (void*)(-mtxfmt),         ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descA2, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, N,   N, 0, 0, N, N, P, Q );
    CHAMELEON_Desc_Create(
        &descQ1, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, M,   N, 0, 0, M, N, P, Q );
    CHAMELEON_Desc_Create(
        &descQ2, CHAMELEON_MAT_ALLOC_TILE, ChamComplexDouble, nb, nb, nb * nb, N,   N, 0, 0, N, N, P, Q );

    CHAMELEON_zplrnt_Tile( descA1, seedA );
    CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., descA2 );

    CHAMELEON_Alloc_Workspace_zgels( M, N, &TS1, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, N, &TT1, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( N, N, &TS2, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( N, N, &TT2, P, Q );

    /*
     * Create the adapted trees to perform the QR factorizations
     */
    {
        libhqr_matrix_t mat = {
            .mt    = descA1->mt,
            .nt    = descA1->nt,
            .nodes = chameleon_desc_datadist_get_iparam(descA1, 0) * chameleon_desc_datadist_get_iparam(descA1, 1),
            .p     = chameleon_desc_datadist_get_iparam(descA1, 0),
        };

        /* Tree for the top matrix */
        libhqr_init_hqr( &qrtreeT, LIBHQR_QR, &mat,
                         -1,        /*low level tree   */
                         -1,        /* high level tree */
                         -1,        /* TS tree size    */
                         chameleon_desc_datadist_get_iparam(descA1, 0), /* High level size */
                         -1,        /* Domino */
                         0          /* TSRR (unstable) */ );

        /* Tree for the bottom matrix */
        mat.mt = descA2->mt;
        mat.nt = descA2->nt;
        libhqr_init_tphqr(
            &qrtreeB, LIBHQR_TSQR, mat.mt, zqdwh_opt_id ? (mat.nt - 1) : 0, &mat,
            /* high level tree (Could be greedy, but flat should reduce the volume of comm) */
            LIBHQR_FLAT_TREE,
            -1,       /* TS tree size    */
            chameleon_desc_datadist_get_iparam(descA2, 0) /* High level size */ );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        assert( 0 );
    }
    else {
        hres = CHAMELEON_zgepdf_qr_Tile( 1, 1, &qrtreeT, &qrtreeB,
                                         descA1, TS1, TT1, descQ1,
                                         descA2, TS2, TT2, descQ2 );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgepdf_qr( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    CHAMELEON_Dealloc_Workspace( &TS1 );
    CHAMELEON_Dealloc_Workspace( &TS2 );
    CHAMELEON_Dealloc_Workspace( &TT1 );
    CHAMELEON_Dealloc_Workspace( &TT2 );

    libhqr_finalize( &qrtreeT );
    libhqr_finalize( &qrtreeB );

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA01, *descA02;
        descA01 = CHAMELEON_Desc_Copy( descA1, CHAMELEON_MAT_ALLOC_TILE );
        descA02 = descA2; /* A2 is useless now */

        CHAMELEON_zplrnt_Tile( descA01, seedA );
        CHAMELEON_zlaset_Tile( ChamUpperLower, 0., 1., descA02 );

        hres += check_zgepdf_qr( args, descA01, descA02, descQ1, descQ2, descA1 );

        CHAMELEON_Desc_Destroy( &descA01 );
    }

    CHAMELEON_Desc_Destroy( &descA1 );
    CHAMELEON_Desc_Destroy( &descA2 );
    CHAMELEON_Desc_Destroy( &descQ1 );
    CHAMELEON_Desc_Destroy( &descQ2 );

    return hres;
}

testing_t   test_zgepdf_qr;
const char *zgepdf_qr_params[] = { "mtxfmt", "nb", "ib", "m", "n", "lda", "seedA", NULL };
const char *zgepdf_qr_output[] = { NULL };
const char *zgepdf_qr_outchk[] = { "||A||", "||A-fact(A)||", "||I-QQ'||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgepdf_qr_init( void ) __attribute__( ( constructor ) );
void
testing_zgepdf_qr_init( void )
{
    test_zgepdf_qr.name   = "zgepdf_qr";
    test_zgepdf_qr.helper = "Polar decomposition factorization with QDWH algorithm";
    test_zgepdf_qr.params = zgepdf_qr_params;
    test_zgepdf_qr.output = zgepdf_qr_output;
    test_zgepdf_qr.outchk = zgepdf_qr_outchk;
    test_zgepdf_qr.fptr_desc = testing_zgepdf_qr_desc;
    test_zgepdf_qr.fptr_std  = NULL;
    test_zgepdf_qr.next   = NULL;

    testing_register( &test_zgepdf_qr );
}
