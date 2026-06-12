/**
 *
 * @file pzhered.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhered parallel algorithm
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Ana Hourcau
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> z d
 *
 */
// ALLOC_WS :  A->mb
// ALLOC_WS :  A->nb
// WS_ADD :  A->mb + A->nb
#include "control/common.h"
#include <coreblas/lapacke.h>

#define A(m, n) A, (m), (n)
#define W(desc, m, n) (desc), (m), (n)

static inline void
chameleon_pzhered_frb( cham_trans_t      trans,
                       cham_uplo_t       uplo,
                       CHAM_desc_t      *A,
                       CHAM_desc_t      *Wnorm,
                       CHAM_desc_t      *Welt,
                       RUNTIME_option_t *options )
{
    double alpha = 1.0;
    double beta = 0.0;

    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = chameleon_desc_datadist_get_iparam( A, 0 );
    int Q  = chameleon_desc_datadist_get_iparam( A, 1 );

    /* Initialize workspaces for tile norms */
    for (m = 0; m < Wnorm->mt; m++)
    {
        for (n = 0; n < NT; n++)
        {
            INSERT_TASK_dlaset(
                options,
                ChamUpperLower, Wnorm->mb, Wnorm->nb,
                alpha, beta,
                W(Wnorm, m, n));
        }
    }

    /* Initialize workspaces */
    for (m = 0; m < Welt->mt; m++)
    {
        for (n = 0; n < Welt->nt; n++)
        {
            INSERT_TASK_dlaset(
                options,
                ChamUpperLower, Welt->mb, Welt->nb,
                alpha, beta,
                W(Welt, m, n));
        }
    }

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for (m = 0; m < MT; m++) {
        int nmin = (uplo == ChamUpper) ? m : 0;
        int nmax = (uplo == ChamLower) ? chameleon_min(m + 1, NT) : NT;

        int tempmm = A->get_blkdim( A, m, DIM_m, M );

        for (n = nmin; n < nmax; n++) {
            int tempnn = A->get_blkdim( A, n, DIM_n, N );

            if ( n == m ) {
                if ( trans == ChamConjTrans ) {
                    INSERT_TASK_zhessq(
                        options, ChamEltwise, uplo, tempmm,
                        A(m, n), W( Wnorm, m, n) );
                }
                else {
                    INSERT_TASK_zsyssq(
                        options, ChamEltwise, uplo, tempmm,
                        A(m, n), W( Wnorm, m, n) );
                }
            }
            else {
                INSERT_TASK_zgessq(
                    options, ChamEltwise, tempmm, tempnn,
                    A(m, n), W( Wnorm, m, n ) );
                INSERT_TASK_zgessq(
                    options, ChamEltwise, tempmm, tempnn,
                    A(m, n), W( Wnorm, n, m ) );
            }
        }
    }

    for(m = 0; m < MT; m++) {
        for(n = Q; n < NT; n++) {
            INSERT_TASK_dplssq(
                options, ChamEltwise, 1, 1, W( Wnorm, m, n), W( Welt, m, n%Q) );
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt( m, 0..Q-1 ) )
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

    /* Compute the norm of each tile, and the full norm */
    for (m = 0; m < MT; m++) {
        int nmin = (uplo == ChamUpper) ? m : 0;
        int nmax = (uplo == ChamLower) ? chameleon_min(m + 1, NT) : NT;

        for (n = nmin; n < nmax; n++) {
            /* Compute the final norm of the tile */
            INSERT_TASK_dplssq2(
                options, 1, W( Wnorm, m, n ) );
        }
    }
    INSERT_TASK_dplssq2(
        options, 1, W( Welt, 0, 0) );

    /**
     * Broadcast the result
     */
    for (m = 0; m < P; m++) {
        for (n = 0; n < Q; n++) {
            if ( ( m != 0 ) || ( n != 0 ) ) {
                INSERT_TASK_dlacpy(
                    options,
                    ChamUpperLower, 1, 1,
                    W( Welt, 0, 0 ), W( Welt, m, n ) );
            }
        }
    }
}

/**
 *
 */
void chameleon_pzhered( cham_trans_t        trans,
                        cham_uplo_t         uplo,
                        double              prec,
                        CHAM_desc_t        *A,
                        RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
    double gnorm, threshold, eps, eps_diag, threshold_diag;

    int workmt, worknt;
    int m, n, P, Q;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    P = chameleon_desc_datadist_get_iparam( A, 0 );
    Q = chameleon_desc_datadist_get_iparam( A, 1 );
    workmt = chameleon_max( A->mt, P );
    worknt = chameleon_max( A->nt, Q );

    RUNTIME_options_ws_alloc( &options, 1, 0 );

    /* Matrix to store the norm of each element */
    chameleon_desc_init( chamctxt, &Wcol, "HERED_Wcol", CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble,
                         2, 1, A->mt * 2, A->nt, A->mt * 2, A->nt, P, Q,
                         NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );

    /* Matrix to compute the global frobenius norm */
    chameleon_desc_init_2dlap( chamctxt, &Welt, "HERED_Welt", ChamRealDouble,
                               2, 1, workmt*2, worknt, P, Q );

    chameleon_pzhered_frb( trans, uplo, A, &Wcol, &Welt, &options );
    Welt.sync = 1; /* Need to be sync to read the norm */

    CHAMELEON_Desc_Flush( &Wcol, sequence );
    CHAMELEON_Desc_Flush( &Welt, sequence );
    CHAMELEON_Desc_Flush( A,     sequence );

    RUNTIME_sequence_wait( chamctxt, sequence );

    gnorm = *((double *)Welt.get_blkaddr( &Welt, A->myrank / Q, A->myrank % Q ));
    chameleon_desc_destroy( &Welt );

    /**
     * Reduce the precision of the tiles if possible
     */
    eps_diag = CHAMELEON_slamch();
    if ( prec < 0. ) {
        eps = CHAMELEON_dlamch();
    }
    else {
        eps = prec;
    }
    threshold = (eps * gnorm) / (double)(chameleon_min(A->mt, A->nt));
    threshold_diag = ( eps < eps_diag ) ? threshold : (eps_diag * gnorm) / (double)(chameleon_min(A->mt, A->nt));

#if defined(CHAMELEON_DEBUG_GERED)
    fprintf( stderr,
             "[%2d] The norm of A is:           %e\n"
             "[%2d] The requested precision is: %e\n"
             "[%2d] The computed threshold is:  %e\n"
             "[%2d] The threshold diag is:      %e\n",
             A->myrank, gnorm,
             A->myrank, eps,
             A->myrank, threshold,
             A->myrank, threshold_diag );
#endif

    for (m = 0; m < A->mt; m++)
    {
        int tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        int nmin   = ( uplo == ChamUpper ) ? m                         : 0;
        int nmax   = ( uplo == ChamLower ) ? chameleon_min(m+1, A->nt) : A->nt;

        for (n = nmin; n < nmax; n++)
        {
            int tempnn = A->get_blkdim( A, n, DIM_n, A->n );

            /*
             * u_{high} = 1e-16 (later should be application accuracy)
             * u_{low} = 1e-8
             * ||A_{i,j}||_F  < u_{high} * || A ||_F / (nt * u_{low})
             * ||A_{i,j}||_F  < threshold / u_{low}
             */
            if ( m == n ) {
                INSERT_TASK_zgered( &options, threshold_diag,
                                    tempmm, tempnn, A( m, n ), W( &Wcol, m, n ) );
            }
            else {
                INSERT_TASK_zgered( &options, threshold,
                                    tempmm, tempnn, A( m, n ), W( &Wcol, m, n ) );
            }
        }
    }
    A->sync = 1; /* A has been updated */

    CHAMELEON_Desc_Flush( A, sequence );
    RUNTIME_sequence_wait( chamctxt, sequence );

    chameleon_desc_destroy( &Wcol );
    RUNTIME_options_ws_free( &options );
    RUNTIME_options_finalize( &options, chamctxt );
}
