/**
 *
 * @file testing_zprint.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zprint testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Lionel Eyraud-Dubois
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
testing_zprint_desc( run_arg_list_t *args, int check )
{
    int rc, hres = 0;

    /* Read arguments */
    intptr_t   mtxfmt    = parameters_getvalue_int( "mtxfmt" );
    int        nb        = run_arg_get_nb( args );
    int        P         = parameters_getvalue_int( "P" );
    int        N         = run_arg_get_int( args, "N", 1000 );
    int        M         = run_arg_get_int( args, "M", N );
    int        LDA       = run_arg_get_int( args, "LDA", M );
    int        l1        = run_arg_get_int( args, "l1", nb / 2 );
    int        l2        = run_arg_get_int( args, "l2", l1 / 3 );
    int        l3        = run_arg_get_int( args, "l3", l2 / 2 );
    cham_rec_t rec       = run_arg_get_rec( args, "rec", ChamRecFull );
    int        rarg      = run_arg_get_int( args, "rarg", 1 );
    int        Q         = parameters_compute_q( P );
    int        list_nb[] = { nb, nb, 0, l2, l3, 0 };
    int        list_mb[] = { M,  nb, 0, l2, l3, 0 };

    /* Descriptors */
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    fprintf( stdout, "--- Tile layout ---\n" );
    rc = CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    if ( rc == CHAMELEON_SUCCESS ) {
        CHAMELEON_Desc_Print( descA );
    }
    else {
        fprintf( stdout, "--- Tile layout (FAILED)---\n" );
        hres++;
    }

    CHAMELEON_Desc_Destroy( &descA );

    fprintf( stdout, "--- Lapack layout ---\n" );
    rc = CHAMELEON_Desc_Create_User(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q,
        chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );

    if ( rc == CHAMELEON_SUCCESS ) {
        CHAMELEON_Desc_Print( descA );
    }
    else {
        fprintf( stdout, "--- Lapack layout (FAILED)---\n" );
        hres++;
    }
    CHAMELEON_Desc_Destroy( &descA );

    if ( CHAMELEON_Comm_size() == 1 )
    {
        fprintf( stdout, "--- Recursive layout (Tile)---\n" );
        rc = CHAMELEON_Recursive_Desc_Create(
            &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
            rec, rarg, list_mb, list_nb, LDA, N, M, N, P, Q,
            NULL, NULL, NULL, NULL, "A" );

        if ( rc == CHAMELEON_SUCCESS ) {
            CHAMELEON_Desc_Print( descA );
        }
        CHAMELEON_Desc_Destroy( &descA );

        fprintf( stdout, "--- Recursive layout (Lapack) ---\n" );
        rc = CHAMELEON_Recursive_Desc_Create(
            &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
            rec, rarg, list_mb, list_nb, LDA, N, M, N, P, Q,
            chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL, "A" );

        if ( rc == CHAMELEON_SUCCESS ) {
            CHAMELEON_Desc_Print( descA );
        }
        CHAMELEON_Desc_Destroy( &descA );
    }

    run_arg_add_fixdbl( args, "time", 1. );
    run_arg_add_fixdbl( args, "gflops", 1. );

    (void)check;
    return hres;
}

testing_t   test_zprint;
const char *zprint_params[] = { "mtxfmt", "nb", "l1", "l2", "l3", "rec", "rarg", "m", "n", "lda", NULL };
const char *zprint_output[] = { NULL };
const char *zprint_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zprint_init( void ) __attribute__( ( constructor ) );
void
testing_zprint_init( void )
{
    test_zprint.name   = "zprint";
    test_zprint.helper = "Print descriptors";
    test_zprint.params = zprint_params;
    test_zprint.output = zprint_output;
    test_zprint.outchk = zprint_outchk;
    test_zprint.fptr_desc = testing_zprint_desc;
    test_zprint.fptr_std  = NULL;
    test_zprint.next   = NULL;

    testing_register( &test_zprint );
}
