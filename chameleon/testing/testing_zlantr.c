/**
 *
 * @file testing_zlantr.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlantr testing
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @author Florent Pruvost
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
#include <coreblas/lapacke.h>
#endif

static cham_fixdbl_t
flops_zlantr( cham_normtype_t ntype, cham_uplo_t uplo, int _M, int _N )
{
    cham_fixdbl_t flops   = 0.;
    cham_fixdbl_t coefabs = 1.;
    cham_fixdbl_t M       = _M;
    cham_fixdbl_t N       = _N;
    cham_fixdbl_t size;
#if defined(PRECISION_z) || defined(PRECISION_c)
    coefabs = 3.;
#endif

    switch ( uplo ) {
        case ChamUpper:
            if ( N > M ) {
                size = ( M * ( M + 1 ) / 2 ) + M * ( N - M );
            }
            else {
                size = N * ( N + 1 ) / 2;
            }
            break;
        case ChamLower:
            if ( M > N ) {
                size = ( N * ( N + 1 ) / 2 ) + N * ( M - N );
            }
            else {
                size = M * ( M + 1 ) / 2;
            }
            break;
        case ChamUpperLower:
        default:
            size = M * N;
    }
    flops = size * coefabs;

    switch ( ntype ) {
        case ChamOneNorm:
            flops += N;
            break;
        case ChamInfNorm:
            flops += M;
            break;
        case ChamMaxNorm:
            break;
        case ChamFrobeniusNorm:
            flops += size;
            break;
        default:;
    }
    (void)flops;
    return size;
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zlantr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int             async     = parameters_getvalue_int( "async" );
    int             nb        = run_arg_get_nb( args );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t     diag      = run_arg_get_diag( args, "diag", ChamNonUnit );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    double       norm;
    CHAM_desc_t *descA;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* Calculates the norm */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlantr_Tile_Async( norm_type, uplo, diag, descA, &norm,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        norm = CHAMELEON_zlantr_Tile( norm_type, uplo, diag, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_znorm( args, ChamTriangular, norm_type, uplo, diag, norm, descA );
    }

    parameters_desc_destroy( &descA );

    return hres;
}
#endif

int
testing_zlantr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    int             api       = parameters_getvalue_int( "api" );
#endif
    int             nb        = run_arg_get_nb( args );
    cham_normtype_t norm_type = run_arg_get_ntype( args, "norm", ChamMaxNorm );
    cham_uplo_t     uplo      = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t     diag      = run_arg_get_diag( args, "diag", ChamNonUnit );
    int             N         = run_arg_get_int( args, "N", 1000 );
    int             M         = run_arg_get_int( args, "M", N );
    int             LDA       = run_arg_get_int( args, "LDA", M );
    int             seedA     = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    double norm;
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrix */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, N, A, LDA, seedA );

    /* Calculates the norm */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    norm = LAPACKE_zlantr( LAPACK_COL_MAJOR, chameleon_lapack_const( norm_type ), uplo, diag, M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        norm = CHAMELEON_zlantr( norm_type, uplo, diag, M, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        norm = CHAMELEON_lapacke_zlantr( CblasColMajor, chameleon_lapack_const( norm_type ), chameleon_lapack_const(uplo), chameleon_lapack_const(diag), M, N, A, LDA );
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
    testing_stop( &test_data, flops_zlantr( norm_type, uplo, M, N ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the solution */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres = check_znorm_std( args, ChamTriangular, norm_type, uplo, diag, norm, M, N, A, LDA );
    }
#endif

    free( A );

    (void)norm;
    (void)check;
    return hres;
}

testing_t   test_zlantr;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zlantr_params[] = { "norm", "uplo",  "diag",
                                "m",      "n",  "lda",  "seedA", NULL };
#else
const char *zlantr_params[] = { "mtxfmt", "nb", "norm", "uplo",  "diag",
                                "m",      "n",  "lda",  "seedA", NULL };
#endif
const char *zlantr_output[] = { NULL };
const char *zlantr_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlantr_init( void ) __attribute__( ( constructor ) );
void
testing_zlantr_init( void )
{
    test_zlantr.name   = "zlantr";
    test_zlantr.helper = "Triangular matrix norm";
    test_zlantr.params = zlantr_params;
    test_zlantr.output = zlantr_output;
    test_zlantr.outchk = zlantr_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zlantr.fptr_desc = NULL;
#else
    test_zlantr.fptr_desc = testing_zlantr_desc;
#endif
    test_zlantr.fptr_std  = testing_zlantr_std;
    test_zlantr.next   = NULL;

    testing_register( &test_zlantr );
}
