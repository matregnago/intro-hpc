/**
 *
 * @file testing_zgesvd.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesvd testing
 *
 * @version 1.4.0
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

static cham_fixdbl_t
flops_zgesvd( int M, int N, int K, cham_job_t jobu, cham_job_t jobvt )
{
    cham_fixdbl_t flops = flops_zgebrd( M, N );

    switch ( jobu ) {
        case ChamAllVec:
            flops += flops_zunmqr( ChamLeft, M, N, N );
            break;
        case ChamOVec:
        case ChamSVec:
            flops += flops_zunmqr( ChamLeft, M, K, K );
            break;
        case ChamNoVec:
            break;
        default:
        ;
    }

    switch ( jobvt ) {
        case ChamAllVec:
            flops += flops_zunmlq( ChamRight, M, N, M );
            break;
        case ChamOVec:
        case ChamSVec:
            flops += flops_zunmlq( ChamRight, K, N, K );
            break;
        case ChamNoVec:
            break;
        default:
        ;
    }

    flops *= -1.; /* make it negative as long as the formulae is not complete */
    return flops;
}

int
testing_zgesvd_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int        async  = parameters_getvalue_int( "async" );
    intptr_t   mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int        nb     = run_arg_get_nb( args );
    int        P      = parameters_getvalue_int( "P" );
    int        N      = run_arg_get_int( args, "N", 1000 );
    int        M      = run_arg_get_int( args, "M", N );
    int        K      = chameleon_min( M, N );
    int        LDA    = run_arg_get_int( args, "LDA", M );
    int        seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    double     cond   = run_arg_get_double( args, "cond", 1.e16 );
    int        mode   = run_arg_get_int( args, "mode", 4 );
    int        Q      = parameters_compute_q( P );
    cham_job_t jobu   = run_arg_get_job( args, "jobu", ChamNoVec );
    cham_job_t jobvt  = run_arg_get_job( args, "jobvt", ChamNoVec );
    int        runtime;

    /* Descriptors */
    CHAM_desc_t           *descA, *descT, *descA0;
    CHAMELEON_Complex64_t *U, *Vt = NULL;
    double                *S, *D;
    int                    LDU   = M;
    int                    LDVt  = N;
    int                    Un;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    CHAMELEON_Get( CHAMELEON_RUNTIME, &runtime );
    if ( runtime == RUNTIME_SCHED_PARSEC ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The SVD is not supported with PaRSEC\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgesvd( M, N, &descT, P, Q );

    if ( (jobu == ChamAllVec) || (jobu == ChamSVec) ) {
        Un  = ( jobu == ChamSVec ) ? K : M;
        U   = malloc( sizeof(CHAMELEON_Complex64_t) * LDU*Un );
    }
    else {
        U = NULL;
    }

    if ( (jobvt == ChamAllVec) || (jobvt == ChamSVec) ) {
        LDVt = ( jobvt == ChamSVec ) ? K : N;
        Vt   = malloc( sizeof(CHAMELEON_Complex64_t) * LDVt*N );
    }
    else {
        Vt = NULL;
    }

    /* Generate the diagonal of eigen/singular values */
    D = malloc( sizeof(double) * K );
    S = malloc( sizeof(double) * K );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms_Tile( ChamDistUniform, seedA, ChamNonsymPosv, D, mode, cond, 1., descA );
    if ( hres != 0 ) {
        free( D );
        free( S );
        free( U );
        free( Vt );
        return hres;
    }
    /*
     * descA0 is defined here because of the cost of zlatms. To copy descA in descA0
     * now prevents to call it again later in the check (indeed descA is modified
     * with the call to CHAMELEON_zgepdf_qdwh_Tile[_Async]).
     */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zlacpy_Tile( ChamUpperLower, descA, descA0 );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgesvd_Tile_Async( jobu, jobvt, descA, S, descT, U, LDU, Vt, LDVt, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descT, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgesvd_Tile( jobu, jobvt, descA, S, descT, U, LDU, Vt, LDVt );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgesvd( M, N, K, jobu, jobvt ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        hres += check_zgesvd( args, jobu, jobvt, descA0, descA, D, S, U, LDU, Vt, LDVt );
        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descT );
    free( S );
    free( D );
    free( U );
    free( Vt );

    return hres;
}

int
testing_zgesvd_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = TESTINGS_TESTDATA_INITIALIZER;
    int        hres      = 0;

    /* Read arguments */
    int        nb    = run_arg_get_nb( args );
    int        N     = run_arg_get_int( args, "N", 1000 );
    int        M     = run_arg_get_int( args, "M", N );
    int        K     = chameleon_min( M, N );
    int        LDA   = run_arg_get_int( args, "LDA", M );
    int        seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    double     cond  = run_arg_get_double( args, "cond", 1.e16 );
    int        mode  = run_arg_get_int( args, "mode", 4 );
    cham_job_t jobu  = run_arg_get_job( args, "jobu", ChamNoVec );
    cham_job_t jobvt = run_arg_get_job( args, "jobvt", ChamNoVec );
    int        runtime;

    /* Descriptors */
    CHAM_desc_t           *descT;
    CHAMELEON_Complex64_t *A, *A0, *U, *Vt;
    double                *S, *D;
    int                    LDU   = M;
    int                    LDVt  = N;
    int                    Un;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    CHAMELEON_Get( CHAMELEON_RUNTIME, &runtime );
    if ( runtime == RUNTIME_SCHED_PARSEC ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: The SVD is not supported with PaRSEC\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    A = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
    CHAMELEON_Alloc_Workspace_zgesvd( M, N, &descT, 1, 1 );

    if ( (jobu == ChamAllVec) || (jobu == ChamSVec) ) {
        Un  = ( jobu == ChamSVec ) ? K : M;
        U   = malloc( sizeof(CHAMELEON_Complex64_t) * LDU*Un );
    }
    else {
        U = NULL;
    }

    if ( (jobvt == ChamAllVec) || (jobvt == ChamSVec) ) {
        LDVt = ( jobvt == ChamSVec ) ? K : N;
        Vt   = malloc( sizeof(CHAMELEON_Complex64_t) * LDVt*N );
    }
    else {
        Vt = NULL;
    }

    /* Generate the diagonal of eigen/singular values */
    D = malloc( sizeof(double) * K );
    S = malloc( sizeof(double) * K );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms( M, N, ChamDistUniform, seedA, ChamNonsymPosv, D, mode, cond, 1., A, LDA );
    if ( hres != 0 ) {
        free( D );
        free( S );
        free( U );
        free( Vt );
        return hres;
    }
    /*
     * A0 is defined here because of the cost of zlatms. To copy A in A0
     * now prevents to call it again later in the check (indeed A is modified
     * with the call to CHAMELEON_zgepdf_qdwh).
     */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {
        A0 = malloc( sizeof(CHAMELEON_Complex64_t) * LDA*N );
        CHAMELEON_zlacpy( ChamUpperLower, M, N, A, LDA, A0, LDA );
    }

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgesvd( jobu, jobvt, M, N, A, LDA, S, descT, U, LDU, Vt, LDVt );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgesvd( M, N, K, jobu, jobvt ) );
    hres = ( hres == CHAMELEON_SUCCESS ) ? 0 : 1;

    /* Checks the factorisation and residue */
    if ( ( hres == CHAMELEON_SUCCESS ) && check ) {

        hres += check_zgesvd_std( args, jobu, jobvt, M, N, A0, A, LDA, D, S, U, LDU, Vt, LDVt );
        free( A0 );
    }

    free( A );
    free( D );
    free( S );
    free( U );
    free( Vt );
    CHAMELEON_Desc_Destroy( &descT );

    return hres;
}

testing_t   test_zgesvd;
const char *zgesvd_params[] = { "mtxfmt", "nb", "jobu", "jobvt", "m", "n", "lda", "seedA", NULL };
const char *zgesvd_output[] = { NULL };
const char *zgesvd_outchk[] = { "||R||", "||I-QQ'||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgesvd_init( void ) __attribute__( ( constructor ) );
void
testing_zgesvd_init( void )
{
    test_zgesvd.name   = "zgesvd";
    test_zgesvd.helper = "Singular Value Decomposition";
    test_zgesvd.params = zgesvd_params;
    test_zgesvd.output = zgesvd_output;
    test_zgesvd.outchk = zgesvd_outchk;
    test_zgesvd.fptr_desc = testing_zgesvd_desc;
    test_zgesvd.fptr_std  = testing_zgesvd_std;
    test_zgesvd.next   = NULL;

    testing_register( &test_zgesvd );
}
