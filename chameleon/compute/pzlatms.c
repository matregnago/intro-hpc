/**
 *
 * @file pzlatms.c
 *
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatms parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/random.h>
#include <coreblas.h>
#endif

#define A(m, n) A,  m,  n

static inline int
zlaset_diag_cpu( void *op_args,
                 cham_uplo_t uplo, int m, int n, int ndata,
                 const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    const double          *D = (const double *)op_args;
    CHAMELEON_Complex64_t *A = CHAM_tile_get_ptr( tileA );

    int tempmm = descA->get_blkdim( descA, m, DIM_m, descA->m );
    int tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );
    int minmn  = chameleon_min( tempmm, tempnn );
    int lda    = tileA->ld;
    int i;

    if ( ndata > 1 ) {
        fprintf( stderr, "zlaset_diag_cpu: supports only one piece of data and %d have been given\n", ndata );
    }

    assert( m == n );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );

    /* Shift to the values corresponding to the tile */
    D += m * descA->mb;

    for( i=0; i<minmn; i++, D++ ) {
        *A = *D;
        A += lda + 1;
    }

    (void)uplo;
    return 0;
}

static cham_map_operator_t zlaset_diag_map = {
    .name        = "zlaset_diag",
    .cpufunc     = zlaset_diag_cpu,
    .cudafunc    = NULL,
    .hipfunc     = NULL,
    .synchronous = 0,
};

/**
 *  Parallel scale of a matrix A
 */
void chameleon_pzlatms( cham_dist_t idist, unsigned long long int seed, cham_sym_t sym,
                        double *D, int mode, double cond, double dmax, CHAM_desc_t *A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t     descU, descV;
    CHAM_desc_t    *descDptr = NULL;
    CHAM_desc_t     descTS, descTT, descD;
    libhqr_matrix_t mat;
    libhqr_tree_t   qrtree;
    cham_trans_t    trans = ChamConjTrans;

    int n, ib, alloc_d = 0;
    int kt    = chameleon_min( A->mt, A->nt );
    int minmn = chameleon_min( A->m,  A->n  );

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    ib = CHAMELEON_IB;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /* Start initialiazing A */
    chameleon_pzlaset( ChamUpperLower, 0., 0., A, sequence, request );

    /* Check if we need to perform UDV or UDU' */
    if ( sym == ChamNonsymPosv ) {
        trans = ChamNoTrans;
    }
    else if ( sym == ChamSymPosv ) {
        trans = ChamTrans;
    }
    else if ( sym == ChamHermPoev ) {
        trans = ChamConjTrans;
    }

#if !defined(CHAMELEON_SIMULATION)
    /* Compute the diagonal D */
    {
        int irsign, rc;

        /* Check what to do with the sign in dlatm1 */
#if defined(PRECISION_z) || defined(PRECISION_c)
        irsign = ( sym == ChamHermGeev );
#else
        irsign = ( sym == ChamHermGeev ) || ( sym == ChamSymPosv );
#endif

        if ( D == NULL ) {
            D = malloc( sizeof(double) * minmn );
            alloc_d = 1;
        }
        rc = CORE_dlatm1( mode, cond, irsign, idist, seed, D, minmn );
        if ( rc != CHAMELEON_SUCCESS ) {
            chameleon_error( "CHAMELEON_zlatms", "Could not generated D. dlatm1 failed" );
            CHAMELEON_Desc_Flush( A, sequence );
            CHAMELEON_Sequence_Wait( sequence );

            if ( alloc_d ) {
                free( D );
            }
            return;
        }

        /* Shift the seed for future use to generate the random unitary matrices */
        seed = CORE_rnd64_jump( 2 * minmn, seed );

        /* Scale by dmax */
        if ( ( mode != 0 ) && ( abs(mode) != 6 ) )
        {
            int    imax = cblas_idamax( minmn, D, 1 );
            double temp = D[imax];

            if ( temp > 0. ) {
                cblas_dscal( minmn, dmax / temp, D, 1 );
            }
            else {
                chameleon_error( "CHAMELEON_zlatms", "Could not scale D (max sing. value is 0.)" );
                CHAMELEON_Desc_Flush( A, sequence );
                CHAMELEON_Sequence_Wait( sequence );

                if ( alloc_d ) {
                    free( D );
                }
                return;
            }
        }
    }
#endif

    /* Copy D to the diagonal of A */
    {
        cham_map_data_t data = {
            .access = ChamRW,
            .desc   = A,
        };
        for (n = 0; n < kt; n++) {
            INSERT_TASK_map(
                &options, ChamUpperLower, n, n,
                1, &data, &zlaset_diag_map, D );
        }
    }

    /**
     * Generate A = U D V
     *
     * U and V are random unitary matrices
     * V = U^t or U^h if A is respectively symmetric or hermitian
     *
     */

    /* Generate U and apply it */
    {
        /* U is of size A->m by min(A->m, A->n) */
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &descU, A->m, minmn );

        chameleon_pzplrnt( &descU, descU.m, 0, 0, seed, sequence, request );

        /* Shift the seed to generate the next random unitary matrix */
#if !defined(CHAMELEON_SIMULATION)
        seed = CORE_rnd64_jump( 2 * minmn * A->m, seed );
#endif

        /* Apply a QR factorization */
        mat.mt    = descU.mt;
        mat.nt    = descU.nt;
        mat.nodes = chameleon_desc_datadist_get_iparam(&descU, 0) * chameleon_desc_datadist_get_iparam(&descU, 1);
        mat.p     = chameleon_desc_datadist_get_iparam(&descU, 0);

        libhqr_init_hqr( &qrtree, LIBHQR_QR, &mat,
                         -1, /*low level tree   */
                         -1, /* high level tree */
                         -1, /* TS tree size    */
                         mat.p,  /* High level size */
                         -1, /* Domino */
                         0   /* TSRR (unstable) */ );

#if defined(CHAMELEON_COPY_DIAG)
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &descD, A->m, minmn );
        descDptr = &descD;
#endif

        chameleon_desc_init_2dtile( chamctxt, &descTS, "LATMS_U_TS", ChamComplexDouble,
                                    ib, descU.nb, ib * descU.mt, descU.nb * descU.nt,
                                    chameleon_desc_datadist_get_iparam(&descU, 0),
                                    chameleon_desc_datadist_get_iparam(&descU, 1) );
        chameleon_desc_init_2dtile( chamctxt, &descTT, "LATMS_U_TT", ChamComplexDouble,
                                    ib, descU.nb, ib * descU.mt, descU.nb * descU.nt,
                                    chameleon_desc_datadist_get_iparam(&descU, 0),
                                    chameleon_desc_datadist_get_iparam(&descU, 1) );

        /* U <= qr(U) */
        chameleon_pzgeqrf_param( 1, kt, &qrtree, &descU,
                                 &descTS, &descTT, descDptr, sequence, request );

        /* A <= U * D */
        chameleon_pzunmqr_param( 0, &qrtree, ChamLeft, ChamNoTrans,
                                 &descU, A, &descTS, &descTT, descDptr, sequence, request );

        if ( trans != ChamNoTrans ) {
            /* A <= (U * D) * U' */
            chameleon_pzunmqr_param( 0, &qrtree, ChamRight, trans,
                                     &descU, A, &descTS, &descTT, descDptr, sequence, request );
        }

        CHAMELEON_Desc_Flush( &descU,  sequence );
        CHAMELEON_Desc_Flush( &descTS, sequence );
        CHAMELEON_Desc_Flush( &descTT, sequence );
        if ( descDptr != NULL ) {
            CHAMELEON_Desc_Flush( descDptr, sequence );
        }
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Sequence_Wait( sequence );

        chameleon_desc_destroy( &descU );
        chameleon_desc_destroy( &descTS );
        chameleon_desc_destroy( &descTT );
        if ( descDptr != NULL ) {
            chameleon_desc_destroy( descDptr );
        }

        libhqr_finalize( &qrtree );
    }

    /* Generate V and apply it */
    if ( trans == ChamNoTrans )
    {
        /* V is of size min(A->m, A->n) by A->n */
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &descV, minmn, A->n );

        chameleon_pzplrnt( &descV, descV.m, 0, 0, seed, sequence, request );

        /* Apply a QR factorization */
        mat.mt    = descV.mt;
        mat.nt    = descV.nt;
        mat.nodes = chameleon_desc_datadist_get_iparam(&descV, 0) * chameleon_desc_datadist_get_iparam(&descV, 1);
        mat.p     = chameleon_desc_datadist_get_iparam(&descV, 1);

        libhqr_init_hqr( &qrtree, LIBHQR_LQ, &mat,
                         -1, /*low level tree   */
                         -1, /* high level tree */
                         -1, /* TS tree size    */
                         mat.p,  /* High level size */
                         -1, /* Domino */
                         0   /* TSRR (unstable) */ );

#if defined(CHAMELEON_COPY_DIAG)
        chameleon_zdesc_copy_and_restrict( chamctxt, A, &descD, minmn, A->n );
        descDptr = &descD;
#endif

        chameleon_desc_init_2dtile( chamctxt, &descTS, "LATMS_V_TS", ChamComplexDouble,
                                    ib, descV.nb, ib * descV.mt, descV.nb * descV.nt,
                                    chameleon_desc_datadist_get_iparam(&descV, 0),
                                    chameleon_desc_datadist_get_iparam(&descV, 1) );
        chameleon_desc_init_2dtile( chamctxt, &descTT, "LATMS_V_TT", ChamComplexDouble,
                                    ib, descV.nb, ib * descV.mt, descV.nb * descV.nt,
                                    chameleon_desc_datadist_get_iparam(&descV, 0),
                                    chameleon_desc_datadist_get_iparam(&descV, 1) );

        /* V <= qr(V) */
        chameleon_pzgelqf_param( 1, kt, &qrtree, &descV,
                                 &descTS, &descTT, descDptr, sequence, request );

        /* A <= (U * D) * V */
        chameleon_pzunmlq_param( 0, &qrtree, ChamRight, ChamNoTrans,
                                 &descV, A, &descTS, &descTT, descDptr, sequence, request );

        CHAMELEON_Desc_Flush( &descV,  sequence );
        CHAMELEON_Desc_Flush( &descTS, sequence );
        CHAMELEON_Desc_Flush( &descTT, sequence );
        if ( descDptr != NULL ) {
            CHAMELEON_Desc_Flush( descDptr, sequence );
        }
        CHAMELEON_Desc_Flush( A, sequence );
        CHAMELEON_Sequence_Wait( sequence );

        chameleon_desc_destroy( &descV );
        chameleon_desc_destroy( &descTS );
        chameleon_desc_destroy( &descTT );
        if ( descDptr != NULL ) {
            chameleon_desc_destroy( descDptr );
        }

        libhqr_finalize( &qrtree );
    }
    /* Mark written data for synchronization */
    A->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);

    if ( alloc_d ) {
        free( D );
    }
    (void)descD;
}
