/**
 *
 * @file pzlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
//ALLOC_WS :  A->mb
//ALLOC_WS :  A->nb
//WS_ADD :  A->mb + A->nb
#include "control/common.h"

#define A( m, n )        A,     (m), (n)
#define W( desc, m, n ) (desc), (m), (n)

static inline void
chameleon_pzlange_one( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                       CHAM_desc_t *Wcol, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options )
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = chameleon_desc_datadist_get_iparam(Welt, 0);
    int Q  = chameleon_desc_datadist_get_iparam(Welt, 1);

    /**
     * Step 1:
     *  For j in [1,P], W(i, n) = reduce( A(i+k*P, n) )
     */
    for(n = 0; n < NT; n++) {
        int mmin = ( uplo == ChamLower ) ? n                      : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, MT) : MT;

        int tempnn = A->get_blkdim( A, n, DIM_n, N );

        for(m = mmin; m < mmax; m++) {
            int tempmm = A->get_blkdim( A, m, DIM_m, M );

            if ( (n == m) && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrasm(
                    options,
                    ChamColumnwise, uplo, diag, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n ) );
            }
            else {
                INSERT_TASK_dzasum(
                    options,
                    ChamColumnwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n ) );
            }

            if ( m >= P ) {
                INSERT_TASK_daxpy(
                    options, tempnn, 1.,
                    W( Wcol, m,   n ), 1,
                    W( Wcol, m%P, n ), 1 );
            }
        }

        /**
         * Step 2:
         *  For each i, W(i, n) = reduce( W(0..P-1, n) )
         */
        for(m = 1; m < P; m++) {
            INSERT_TASK_daxpy(
                options, tempnn, 1.,
                W( Wcol, m, n ), 1,
                W( Wcol, 0, n ), 1 );
        }

        INSERT_TASK_dlange(
            options,
            ChamMaxNorm, 1, tempnn, A->nb,
            W( Wcol, 0, n ),
            W( Welt, 0, n ) );
    }

    /**
     * Step 3:
     *  For n in 0..Q-1, W(m, n) = max( W(m, n..nt[Q] ) )
     */
    for(n = Q; n < NT; n++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, 0, n   ),
            W( Welt, 0, n%Q ) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(n = 1; n < Q; n++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, 0, n ),
            W( Welt, 0, 0 ) );
    }
}

static inline void
chameleon_pzlange_inf( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A,
                       CHAM_desc_t *Wcol, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = chameleon_desc_datadist_get_iparam(Welt, 0);
    int Q  = chameleon_desc_datadist_get_iparam(Welt, 1);

    /**
     * Step 1:
     *  For j in [1,Q], Wcol(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = A->get_blkdim( A, m, DIM_m, M );

        for(n = nmin; n < nmax; n++) {
            int tempnn = A->get_blkdim( A, n, DIM_n, N );

            if ( (n == m)  && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrasm(
                    options,
                    ChamRowwise, uplo, diag, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n) );
            }
            else {
                INSERT_TASK_dzasum(
                    options,
                    ChamRowwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n) );
            }

            if ( n >= Q ) {
                INSERT_TASK_daxpy(
                    options, tempmm, 1.,
                    W( Wcol, m, n   ), 1,
                    W( Wcol, m, n%Q ), 1 );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Wcol(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_daxpy(
                options, tempmm, 1.,
                W( Wcol, m, n ), 1,
                W( Wcol, m, 0 ), 1 );
        }

        INSERT_TASK_dlange(
            options,
            ChamMaxNorm, tempmm, 1, A->nb,
            W( Wcol, m, 0), W( Welt, m, 0));
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Wcol(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, 0, 0) );
    }
}

static inline void
chameleon_pzlange_max( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = chameleon_desc_datadist_get_iparam(Welt, 0);
    int Q  = chameleon_desc_datadist_get_iparam(Welt, 1);

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = A->get_blkdim( A, m, DIM_m, M );

        for(n = nmin; n < nmax; n++) {
            int tempnn = A->get_blkdim( A, n, DIM_n, N );

            if ( (n == m)  && (uplo != ChamUpperLower) ) {
                INSERT_TASK_zlantr(
                    options,
                    ChamMaxNorm, uplo, diag, tempmm, tempnn, A->nb,
                    A(m, n), W( Welt, m, n));
            }
            else {
                INSERT_TASK_zlange(
                    options,
                    ChamMaxNorm, tempmm, tempnn, A->nb,
                    A(m, n), W( Welt, m, n ));
            }

            if ( n >= Q ) {
                INSERT_TASK_dlange_max(
                    options,
                    W( Welt, m, n   ),
                    W( Welt, m, n%Q ) );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dlange_max(
                options,
                W( Welt, m, n ),
                W( Welt, m, 0 ) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m,   0 ),
            W( Welt, m%P, 0 ) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0 ),
            W( Welt, 0, 0 ) );
    }
}

static inline void
chameleon_pzlange_frb( cham_uplo_t uplo, cham_diag_t diag, CHAM_desc_t *A, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = chameleon_desc_datadist_get_iparam(Welt, 0);
    int Q  = chameleon_desc_datadist_get_iparam(Welt, 1);

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = A->get_blkdim( A, m, DIM_m, M );

        for(n = nmin; n < nmax; n++) {
            int tempnn = A->get_blkdim( A, n, DIM_n, N );

            if ( (n == m) && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrssq(
                    options,
                    uplo, diag, tempmm, tempnn,
                    A(m, n), W( Welt, m, n) );
            }
            else {
                INSERT_TASK_zgessq(
                    options,
                    ChamEltwise,
                    tempmm, tempnn,
                    A(m, n), W( Welt, m, n) );
            }

            if ( n >= Q ) {
                INSERT_TASK_dplssq(
                    options, ChamEltwise, 1, 1, W( Welt, m, n), W( Welt, m, n%Q) );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dplssq(
                options, ChamEltwise, 1, 1, W( Welt, m, n), W( Welt, m, 0) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, 1, W( Welt, m, 0), W( Welt, m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, 1, W( Welt, m, 0), W( Welt, 0, 0) );
    }

    INSERT_TASK_dplssq2(
        options, 1, W( Welt, 0, 0) );
}

/**
 *
 */
void chameleon_pzlange_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                                CHAM_desc_t *A, double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
    double alpha = 0.0;
    double beta  = 0.0;

    int workmt, worknt;
    int m, n, wcol_init = 0;

    int P = chameleon_desc_datadist_get_iparam(A, 0);
    int Q = chameleon_desc_datadist_get_iparam(A, 1);

    *result = -1.0;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    workmt = chameleon_max( A->mt, P );
    worknt = chameleon_max( A->nt, Q );

    switch ( norm ) {
    case ChamOneNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        chameleon_desc_init_2dtile( chamctxt, &Wcol, "LANGE_Wcol", ChamRealDouble,
                                    1, A->nb, workmt, worknt * A->nb, P, Q );
        wcol_init = 1;

        /*
         * Use the global allocator for Welt, otherwise flush may free the data before the result is read.
         */
        chameleon_desc_init_2dlap( chamctxt, &Welt, "LANGE_Welt", ChamRealDouble,
                                   1, 1, P, worknt, P, Q );
        break;

        /*
         *  ChamInfNorm
         */
    case ChamInfNorm:
        RUNTIME_options_ws_alloc( &options, A->mb, 0 );

        chameleon_desc_init_2dtile( chamctxt, &Wcol, "LANGE_Wcol", ChamRealDouble,
                                    A->mb, 1, workmt * A->mb, worknt, P, Q );
        wcol_init = 1;

        chameleon_desc_init_2dlap( chamctxt, &Welt, "LANGE_Welt", ChamRealDouble,
                                   1, 1, workmt, Q, P, Q );
        break;

        /*
         *  ChamFrobeniusNorm
         */
    case ChamFrobeniusNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        alpha = 1.;
        chameleon_desc_init_2dlap( chamctxt, &Welt, "LANGE_Welt", ChamRealDouble,
                                   2, 1, workmt*2, worknt, P, Q );
        break;

        /*
         *  ChamMaxNorm
         */
    case ChamMaxNorm:
    default:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        chameleon_desc_init_2dlap( chamctxt, &Welt, "LANGE_Welt", ChamRealDouble,
                                   1, 1, workmt, worknt, P, Q );
    }

    /* Initialize workspaces */
    if ( (norm == ChamInfNorm) ||
         (norm == ChamOneNorm) )
    {
        /* Initialize Wcol tile */
        for(m = 0; m < Wcol.mt; m++) {
            for(n = 0; n < Wcol.nt; n++) {
                INSERT_TASK_dlaset(
                    &options,
                    ChamUpperLower, Wcol.mb, Wcol.nb,
                    alpha, beta,
                    W( &Wcol, m, n ) );
            }
        }
    }
    for(m = 0; m < Welt.mt; m++) {
        for(n = 0; n < Welt.nt; n++) {
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, Welt.mb, Welt.nb,
                alpha, beta,
                W( &Welt, m, n ) );
        }
    }

    switch ( norm ) {
    case ChamOneNorm:
        chameleon_pzlange_one( uplo, diag, A, &Wcol, &Welt, &options );
        CHAMELEON_Desc_Flush( &Wcol, sequence );
        break;

    case ChamInfNorm:
        chameleon_pzlange_inf( uplo, diag, A, &Wcol, &Welt, &options );
        CHAMELEON_Desc_Flush( &Wcol, sequence );
        break;

    case ChamFrobeniusNorm:
        chameleon_pzlange_frb( uplo, diag, A, &Welt, &options );
        break;

    case ChamMaxNorm:
    default:
        chameleon_pzlange_max( uplo, diag, A, &Welt, &options );
    }

    /**
     * Broadcast the result
     */
    for(m = 0; m < chameleon_desc_datadist_get_iparam(A, 0); m++) {
        for(n = 0; n < chameleon_desc_datadist_get_iparam(A, 1); n++) {
            if ( (m != 0) || (n != 0) ) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1,
                    W( &Welt, 0, 0 ), W( &Welt, m, n ) );
            }
        }
    }
    Welt.sync  = 1; /* Need to be sync to read the norm */

    if ( wcol_init ) {
        CHAMELEON_Desc_Flush( &Wcol, sequence );
    }
    CHAMELEON_Desc_Flush( &Welt, sequence );
    CHAMELEON_Desc_Flush( A, sequence );
    RUNTIME_sequence_wait( chamctxt, sequence );

    *result = *((double *)Welt.get_blkaddr( &Welt, A->myrank / chameleon_desc_datadist_get_iparam(A, 1), A->myrank % chameleon_desc_datadist_get_iparam(A, 1) ));

    if ( wcol_init ) {
        chameleon_desc_destroy( &Wcol );
    }
    chameleon_desc_destroy( &Welt );

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
