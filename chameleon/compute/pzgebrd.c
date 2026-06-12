/**
 *
 * @file pzgebrd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgebrd parallel algorithm
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/lapacke.h"
#endif

#define A(m,n) A, (m), (n)

static inline void
chameleon_pzgebrd_ge2gb( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host   = 0;
    int    k, m, n, ib;
    CHAM_desc_t *Abis = NULL;
    CHAM_desc_t *Tbis = NULL;
    CHAM_desc_t *Dbis = NULL;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }
    ib = CHAMELEON_IB;

    /*
     * zgeqrt/zgelqt   = A->nb * (ib+1)
     * zunmqr/zunmlq   = A->nb * ib
     * ztpqrt/ztplqt   = A->nb * (ib+1)
     * ztpmqrt/ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr/zunmlq   =     A->nb * ib
     * ztpmqrt/ztpmlqt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    if ( A->m >= A->n ){
        /*
         * Create submatrix descriptor of the 3 matrices without the first column
         */
        if ( A->nt > 1 ) {
            Abis = chameleon_desc_submatrix( A, 0, A->nb, A->m, chameleon_max( A->n - A->nb, 0 ) );
            Tbis = chameleon_desc_submatrix( T, 0, T->nb, T->m, chameleon_max( T->n - T->nb, 0 ) );
            if ( D != NULL ) {
                Dbis = chameleon_desc_submatrix( D, 0, D->nb, D->m, chameleon_max( D->n - D->nb, 0 ) );
            }
        }

        for ( k = 0; k < A->nt; k++ ) {
            RUNTIME_iteration_push(chamctxt, k);

            chameleon_pzgeqrf_step( genD, k, ib, A, T, D, &options, sequence );

            /* Restore the original location of the tiles */
            for (n = k; n < A->nt; n++) {
                RUNTIME_data_migrate( sequence, A(k, n),
                                      A->get_rankof( A, k, n ) );
            }

            if ( k+1 < A->nt ) {
                chameleon_pzgelqf_step( genD, k, ib, Abis, Tbis, Dbis, &options, sequence );

                /* Restore the original location of the tiles */
                for (m = k; m < A->mt; m++) {
                    RUNTIME_data_migrate( sequence, A( m, k+1 ),
                                          A->get_rankof( A, m, k+1 ) );
                }
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    else {
        /*
         * Create submatrix descriptor of the 3 matrices without the first row
         */
        if ( A->mt > 1 ) {
            Abis = chameleon_desc_submatrix( A, A->mb, 0, chameleon_max( A->m - A->mb, 0 ), A->n );
            Tbis = chameleon_desc_submatrix( T, T->mb, 0, chameleon_max( T->m - T->mb, 0 ), T->n );
            if ( D != NULL ) {
                Dbis = chameleon_desc_submatrix( D, D->mb, 0, chameleon_max( D->m - D->mb, 0 ), D->n );
            }
        }

        for ( k = 0; k < A->mt; k++ ) {
            RUNTIME_iteration_push(chamctxt, k);

            chameleon_pzgelqf_step( genD, k, ib, A, T, D, &options, sequence );

            /* Restore the original location of the tiles */
            for (m = k; m < A->mt; m++) {
                RUNTIME_data_migrate( sequence, A( m, k ),
                                      A->get_rankof( A, m, k ) );
            }

            if ( k+1 < A->mt ) {
                chameleon_pzgeqrf_step( genD, k, ib, Abis, Tbis, Dbis, &options, sequence );

                /* Restore the original location of the tiles */
                for (n = k; n < A->nt; n++) {
                    RUNTIME_data_migrate( sequence, A( k+1, n ),
                                          A->get_rankof( A, k+1, n ) );
                }
            }
            RUNTIME_iteration_pop(chamctxt);
        }
    }

    if ( Abis != NULL ) {
        free( Abis );
    }
    if ( Tbis != NULL ) {
        free( Tbis );
    }
    if ( Dbis != NULL ) {
        free( Dbis );
    }
    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );
    if ( D != NULL ) {
        CHAMELEON_Desc_Flush( D, sequence );
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}

static inline int
chameleon_pzgebrd_gb2bd( cham_job_t jobu, cham_job_t jobvt, CHAM_desc_t *A,
                         CHAMELEON_Complex64_t *U, int LDU,
                         CHAMELEON_Complex64_t *VT, int LDVT,
                         double *E, double *S,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t descAB;
    cham_uplo_t uplo;
    int M, N, MINMN, NB, LDAB, ABn;
#if !defined(CHAMELEON_SIMULATION)
    int info;
    int KL, KU;
    char gbbrd_vect;
#endif

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return sequence->status;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    M     = A->m;
    N     = A->n;
    MINMN = chameleon_min(M, N);
    NB    = A->mb;
    LDAB  = NB + 1;
    uplo  = M >= N ? ChamUpper : ChamLower;
    ABn   = MINMN;

    /* Allocate band structure */
    chameleon_zdesc_alloc( descAB, "GEBRD_AB",
                           LDAB, NB,  /* mb, nb */
                           LDAB, ABn, /* m,  n  */
                           NULL );

    /* Convert matrix to band form */
    chameleon_pztile2band( uplo, A, &descAB, sequence, request );

#if !defined(CHAMELEON_SIMULATION)
    /* NCC = 0, C = NULL, we do not update any matrix with new singular vectors */
    /* On exit, AB = U (S +~ E) VT */
    KL = uplo == ChamUpper ? 0  : NB;
    KU = uplo == ChamUpper ? NB : 0;

    /* Manage the case where only singular values are required */
    if ( jobu == ChamNoVec ) {
        if ( jobvt == ChamNoVec ) {
            gbbrd_vect = 'N';
        }
        else {
            gbbrd_vect = 'P';
        }
    }
    else {
        if ( jobvt == ChamNoVec ) {
            gbbrd_vect = 'Q';
        }
        else {
            gbbrd_vect = 'B';
        }
    }
#endif
    /* AB being modified must be synced */
    descAB.sync = 1;

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( &descAB, sequence );
    chameleon_sequence_wait( chamctxt, sequence );

#if !defined(CHAMELEON_SIMULATION)
    info = LAPACKE_zgbbrd( LAPACK_COL_MAJOR, gbbrd_vect, M, N, 0, KL, KU,
                           (CHAMELEON_Complex64_t *) descAB.mat, LDAB, S, E,
                           U, LDU, VT, LDVT, NULL, 1 );
    if ( info != 0 ) {
        fprintf( stderr, "CHAMELEON_zgesvd_Tile_Async: LAPACKE_zgbbrd = %d\n", info );
    }
    assert( info == 0 );
#endif /* !defined(CHAMELEON_SIMULATION) */

    chameleon_desc_destroy( &descAB );

    RUNTIME_options_finalize( &options, chamctxt );

    return CHAMELEON_SUCCESS;
}

int chameleon_pzgebrd( int genD, cham_job_t jobu, cham_job_t jobvt,
                       CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                       CHAMELEON_Complex64_t *U, int LDU,
                       CHAMELEON_Complex64_t *VT, int LDVT,
                       double *E, double *S,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *subA, *subT, *subUVT, *subD;
    CHAM_desc_t descUl, descUt;
    CHAM_desc_t descVTl, descVTt;
    int M, N, NB;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return sequence->status;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    chameleon_pzgebrd_ge2gb( genD, A, T, D, sequence, request );
    chameleon_pzgebrd_gb2bd( jobu, jobvt, A, U, LDU, VT, LDVT, E, S, sequence, request );

    /* Update U and Vt according to jobu and jobvt */
    subA   = NULL;
    subT   = NULL;
    subUVT = NULL;
    subD   = NULL;
    M      = A->m;
    N      = A->n;
    NB     = A->mb;

    if ( jobu != ChamNoVec ) {
        chameleon_zlap2tile( chamctxt, "U", &descUl, &descUt, ChamDescInout, ChamUpperLower,
                             U, NB, NB, LDU, M, M, M, sequence, request );

        if ( M < N ) {
            subA   = chameleon_desc_submatrix( A,       chameleon_min(A->mb, A->m),             0,
                                                        chameleon_max(0, A->m - A->mb),         A->n );
            subUVT = chameleon_desc_submatrix( &descUt, chameleon_min(descUt.mb, descUt.m),     0,
                                                        chameleon_max(0, descUt.m - descUt.mb), descUt.n);
            subT   = chameleon_desc_submatrix( T,       chameleon_min(T->mb, T->m),             0,
                                                        chameleon_max(0, T->m - T->mb),         T->n );
            if ( D != NULL ) {
                subD = chameleon_desc_submatrix( D,     chameleon_min(D->mb, D->m),             0,
                                                        chameleon_max(0, D->m - D->mb),         D->n );
            }
            chameleon_pzunmqr( 0, ChamLeft, ChamNoTrans, subA, subUVT, subT, subD, sequence, request );

            free( subA );
            free( subUVT );
            free( subT );
            if ( D != NULL ) {
                free( subD );
            }
        }
        else {
            chameleon_pzunmqr( 0, ChamLeft, ChamNoTrans, A, &descUt, T, D, sequence, request );
        }

        chameleon_ztile2lap( chamctxt, &descUl, &descUt, ChamDescInout, ChamUpperLower, sequence, request );
    }
    if ( jobvt != ChamNoVec ) {
        chameleon_zlap2tile( chamctxt, "VT", &descVTl, &descVTt, ChamDescInout, ChamUpperLower,
                             VT, NB, NB, LDVT, N, N, N, sequence, request );

        if ( M < N ){
            chameleon_pzunmlq( 0, ChamRight, ChamNoTrans, A, &descVTt, T, D, sequence, request );
        }
        else {
            subA   = chameleon_desc_submatrix( A,        0,         chameleon_min(A->nb, A->n),
                                                         A->m,      chameleon_max(0, A->n - A->nb) );
            subUVT = chameleon_desc_submatrix( &descVTt, 0,         chameleon_min(descVTt.nb, descVTt.n),
                                                         descVTt.m, chameleon_max(0, descVTt.n - descVTt.nb) );
            subT   = chameleon_desc_submatrix( T,        0,         chameleon_min(T->nb, T->n),
                                                         T->m,      chameleon_max(0, T->n - T->nb) );
            if ( D != NULL ) {
                subD = chameleon_desc_submatrix( D,      0,         chameleon_min(D->nb, D->n),
                                                         D->m,      chameleon_max(0, D->n - D->nb) );
            }

            chameleon_pzunmlq( 0, ChamRight, ChamNoTrans, subA, subUVT, subT, subD, sequence, request );

            free( subA );
            free( subUVT );
            free( subT );
            if ( D != NULL ) {
                free( subD );
            }
        }

        chameleon_ztile2lap( chamctxt, &descVTl, &descVTt,
                             ChamDescInout, ChamUpperLower, sequence, request );
    }
    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( T, sequence );
    if ( D != NULL ) {
        CHAMELEON_Desc_Flush( D, sequence );
    }
    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    if ( jobu != ChamNoVec ) {
        chameleon_ztile2lap_cleanup( chamctxt, &descUl,  &descUt  );
    }
    if ( jobvt != ChamNoVec ) {
        chameleon_ztile2lap_cleanup( chamctxt, &descVTl, &descVTt );
    }

    RUNTIME_options_finalize( &options, chamctxt );

    return CHAMELEON_SUCCESS;
}
