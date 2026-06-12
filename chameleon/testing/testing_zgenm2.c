/**
 *
 * @file testing_zgenm2.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgenm2 testing
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
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#endif

static cham_fixdbl_t
flops_zgenm2( int M, int N )
{
    double coefabs = 1.;
#if defined(PRECISION_z) || defined(PRECISION_c)
    coefabs = 3.;
#endif

    return coefabs * (cham_fixdbl_t)(M) * (cham_fixdbl_t)(N);
}

int
testing_zgenm2_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );
    int      minMN  = chameleon_min( M, N );
    double   cond   = run_arg_get_double( args, "cond", 1.e16 );
    int      mode   = run_arg_get_int( args, "mode", 4 );
    double   tol    = 1.e-1;

    /* Descriptors */
    double       norm;
    CHAM_desc_t *descA;
    double      *D, dmax = 1.;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Generate the diagonal of eigen/singular values */
    D = malloc( sizeof(double) * minMN );
#if !defined(CHAMELEON_SIMULATION)
    hres = CORE_dlatm1( mode, cond, 0, ChamDistUniform, seedA, D, minMN );
    if ( hres != 0 ) {
        free( D );
        return hres;
    }

    /* Save the largest absolute value */
    hres = cblas_idamax( minMN, D, 1 );
    dmax = fabs( D[hres] );
#else
    (void)mode;
#endif

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms_Tile(
        ChamDistUniform, seedA, ChamNonsymPosv, D, 0, cond, 0., descA );
    free( D );
    if ( hres != 0 ) {
        return hres;
    }

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgenm2_Tile_Async( tol, descA, &norm,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        norm = CHAMELEON_zgenm2_Tile( tol, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgenm2( M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        double res = fabs(dmax - norm) / (dmax * tol);

        run_arg_add_double( args, "||A||", dmax );
        run_arg_add_double( args, "||B||", norm );
        run_arg_add_double( args, "||R||", res );

        if ( isnan(res) || isinf(res) || (res > 10.0) ) {
            hres = 1;
        }
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

int
testing_zgenm2_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int    nb    = run_arg_get_nb( args );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int    minMN = chameleon_min( M, N );
    double cond  = run_arg_get_double( args, "cond", 1.e16 );
    int    mode  = run_arg_get_int( args, "mode", 4 );
    double tol   = 1.e-1;

    /* Descriptors */
    double       norm;
    CHAMELEON_Complex64_t *A;
    double      *D, dmax = 1.;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Generate the diagonal of eigen/singular values */
    D = malloc( sizeof(double) * minMN );
#if !defined(CHAMELEON_SIMULATION)
    hres = CORE_dlatm1( mode, cond, 0, ChamDistUniform, seedA, D, minMN );
    if ( hres != 0 ) {
        free( D );
        return hres;
    }

    /* Save the largest absolute value */
    hres = cblas_idamax( minMN, D, 1 );
    dmax = fabs( D[hres] );
#else
    (void)mode;
#endif

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms( M, N, ChamDistUniform, seedA, ChamNonsymPosv, D, 0, cond, 0., A, LDA );
    free( D );
    if ( hres != 0 ) {
        return hres;
    }

    /* Calculates the norm */
    testing_start( &test_data );
    norm = CHAMELEON_zgenm2( tol, M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgenm2( M, N ) );

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        double res = fabs(dmax - norm) / (dmax * tol);

        run_arg_add_double( args, "||A||", dmax );
        run_arg_add_double( args, "||B||", norm );
        run_arg_add_double( args, "||R||", res );

        if ( isnan(res) || isinf(res) || (res > 10.0) ) {
            hres = 1;
        }
    }

    free( A );

    return hres;
}

testing_t   test_zgenm2;
const char *zgenm2_params[] = { "mtxfmt", "nb", "m", "n", "lda", "seedA", "cond", "mode", NULL };
const char *zgenm2_output[] = { NULL };
const char *zgenm2_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgenm2_init( void ) __attribute__( ( constructor ) );
void
testing_zgenm2_init( void )
{
    test_zgenm2.name   = "zgenm2";
    test_zgenm2.helper = "General matrix two-norm estimator";
    test_zgenm2.params = zgenm2_params;
    test_zgenm2.output = zgenm2_output;
    test_zgenm2.outchk = zgenm2_outchk;
    test_zgenm2.fptr_desc = testing_zgenm2_desc;
    test_zgenm2.fptr_std  = testing_zgenm2_std;
    test_zgenm2.next   = NULL;

    testing_register( &test_zgenm2 );
}
