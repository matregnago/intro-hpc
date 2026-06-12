/**
 *
 * @file testing_zlaswp.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp testing
 *
 * @version 1.4.0
 * @author Matteo Marcos
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon/constants.h"
#include "chameleon/struct.h"
#include "testings.h"
#include "chameleon/chameleon_z.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <chameleon/getenv.h>
#include <coreblas/lapacke.h>
#include <chameleon/tasks.h>

static cham_fixdbl_t
flops_zlaswp( cham_fixdbl_t M, cham_fixdbl_t N )
{
    cham_fixdbl_t flops;

    flops = M * N * sizeof( CHAMELEON_Complex64_t );

    return flops;
}

static void
testing_zlaswp_ipiv_gen( int *IPIV,
                         int  k )
{
    int i;

    for ( i = 0; i < k; i++ ) {
        IPIV[i] = testing_ialea() % ( k - i ) + i + 1;
    }
}

int
testing_zlaswp_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int         async   = parameters_getvalue_int( "async" );
    int         nb      = run_arg_get_nb(  args );
    int         P       = parameters_getvalue_int( "P" );
    cham_side_t side    = run_arg_get_side( args, "side", ChamLeft );
    cham_dir_t  dir     = run_arg_get_dir( args,  "dir", ChamDirForward );
    int         N       = run_arg_get_int( args, "N", 1000 );
    int         M       = run_arg_get_int( args, "M", N );
    int         LDA     = run_arg_get_int( args, "LDA", M );
    int         K1      = run_arg_get_int( args, "K1", 1 );
    int         K2      = run_arg_get_int( args, "K2", ( side == ChamLeft ) ? M : N );
    int         seedA   = run_arg_get_int( args, "seedA", testing_ialea() );
    int         Q       = parameters_compute_q( P );

    int  K    = ( side == ChamLeft ) ? M : N;
    int  kb   = nb;
    int *IPIV;

    /* Descriptors */
    CHAM_desc_t *descA;
    CHAM_ipiv_t *descIPIV;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );
    if ( descA->get_rankof_init != chameleon_getrankof_2d ) {
        parameters_desc_destroy( &descA );
        fprintf( stderr,
                 "SKIPPED: zlaswp kernel only supports 2DBC data distributions\n" );
        return -1;
    }

    IPIV = malloc( sizeof(int) * K );
    testing_zlaswp_ipiv_gen( IPIV, K );
    CHAMELEON_Ipiv_Create( &descIPIV, side, kb, K, P, P*Q, IPIV );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* IPIV is initialized with random values that are propagated to the descriptor (should be changed in the future) */

    if ( async ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( side, descA );
    }

    /* Calculates the perumtation */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlaswp_Tile_Async( side, dir, descA, K1, K2, descIPIV, ws,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA,    test_data.sequence );
        CHAMELEON_Ipiv_Flush( descIPIV, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlaswp_Tile( side, dir, descA, K1, K2, descIPIV );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlaswp( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    if ( ws != NULL ) {
        CHAMELEON_zlaswp_WS_Free( ws );
    }

#if !defined(CHAMELEON_SIMULATION)
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        CHAM_desc_t *descA0, *descA0c;
        int          INCX = ( dir == ChamDirForward ) ? 1 : -1;

        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_Desc_Create_User(
            &descA0c, (void*)CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
            nb, nb, nb*nb, M, N, 0, 0, M, N, 1, 1,
            chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );

        CHAMELEON_zplrnt_Tile( descA0c, seedA );

        if ( CHAMELEON_Comm_rank() == 0 ) {
            int nbelt   = ( side == ChamLeft ) ? N : M;
            int storage = ( side == ChamLeft ) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
            LAPACKE_zlaswp( storage, nbelt, descA0c->mat, M, K1, K2, IPIV, INCX );
        }

        CHAMELEON_zlacpy_Tile( ChamUpperLower, descA0c, descA0 );
        CHAMELEON_Desc_Destroy( &descA0c );

        hres += check_zmatrices( args, ChamUpperLower, descA, descA0 );

        CHAMELEON_Desc_Destroy( &descA0 );
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    CHAMELEON_Ipiv_Destroy( &descIPIV );
    parameters_desc_destroy( &descA );
    free( IPIV );

    return hres;
}

testing_t   test_zlaswp;
const char *zlaswp_params[] = { "mtxfmt", "nb", "side", "dir", "m", "n", "lda", "k1", "k2", "seedA", NULL };
const char *zlaswp_output[] = { NULL };
const char *zlaswp_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlaswp_init( void ) __attribute__( ( constructor ) );
void
testing_zlaswp_init( void )
{
    test_zlaswp.name      = "zlaswp";
    test_zlaswp.helper    = "Row interchange on general matrices";
    test_zlaswp.params    = zlaswp_params;
    test_zlaswp.output    = zlaswp_output;
    test_zlaswp.outchk    = zlaswp_outchk;
    test_zlaswp.fptr_desc = testing_zlaswp_desc;
    test_zlaswp.next      = NULL;

    testing_register( &test_zlaswp );
}

