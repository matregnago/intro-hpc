/**
 *
 * @file pzhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(  _m_, _n_ ) A,  (_m_), (_n_)
#define B(  _m_, _n_ ) B,  (_m_), (_n_)
#define C(  _m_, _n_ ) C,  (_m_), (_n_)
#define WA( _m_, _n_ ) WA, (_m_), (_n_)
#define WB( _m_, _n_ ) WB, (_m_), (_n_)

/**
 *  Parallel tile matrix-matrix multiplication.
 *  Generic algorithm for any data distribution with a stationnary A.
 *
 * Assuming A has been setup with a proper getrank function to account for symmetry
 */
static inline void
chameleon_pzhemm_Astat( CHAM_context_t *chamctxt, cham_side_t side, cham_uplo_t uplo,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        RUNTIME_option_t *options )
{
    const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;
    RUNTIME_sequence_t *sequence = options->sequence;
    RUNTIME_request_t  *request  = options->request;
    int                 k, m, n, l, Am, An;
    int                 tempmm, tempnn, tempkn, tempkm;
    int                 myrank = RUNTIME_comm_rank( chamctxt );
    int8_t             *reduceC = calloc( C->mt * C->nt, sizeof(int8_t) );

    /* Set C tiles to redux mode */
    for (n = 0; n < C->nt; n++) {
        for (m = 0; m < C->mt; m++) {
            /* The node owns the C tile. */
            if ( C->get_rankof( C(m, n) ) == myrank ) {
                reduceC[ n * C->mt + m ] = 1;
                RUNTIME_zgersum_set_methods( C(m, n) );
                continue;
            }

            /*
             * The node owns the A tile that will define the locality of the
             * computations.
             */
            /* Select row or column based on side */
            l = ( side == ChamLeft ) ? m : n;

            if ( uplo == ChamLower ) {
                for (k = 0; k < A->mt; k++) {
                    Am = k;
                    An = k;

                    if (k < l) {
                        Am = l;
                    }
                    else if (k > l) {
                        An = l;
                    }

                    if ( A->get_rankof( A( Am, An ) ) == myrank ) {
                        reduceC[ n * C->mt + m ] = 1;
                        RUNTIME_zgersum_set_methods( C(m, n) );
                        break;
                    }
                }
            }
            else {
                for (k = 0; k < A->mt; k++) {
                    Am = k;
                    An = k;

                    if (k < l) {
                        An = l;
                    }
                    else if (k > l) {
                        Am = l;
                    }

                    if ( A->get_rankof( A( Am, An ) ) == myrank ) {
                        reduceC[ n * C->mt + m ] = 1;
                        RUNTIME_zgersum_set_methods( C(m, n) );
                        break;
                    }
                }
            }
        }
    }

    for(n = 0; n < C->nt; n++) {
        tempnn = C->get_blkdim( C, n, DIM_n, C->n );
        for(m = 0; m < C->mt; m++) {
            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

            /* Scale C */
            options->forcesub = 0;
            INSERT_TASK_zlascal( options, ChamUpperLower, tempmm, tempnn, C->mb,
                                 beta, C, m, n );
            options->forcesub = reduceC[ n * C->mt + m ];

            /*
             *  ChamLeft / ChamLower
             */
            if (side == ChamLeft) {
                if (uplo == ChamLower) {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = C->get_blkdim( C, k, DIM_m, C->m );

                        if (k < m) {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k),  /* lda * K */
                                       B(k, n),  /* ldb * Y */
                                zone,  C(m, n)); /* ldc * Y */
                        }
                        else if (k == m) {
                                INSERT_TASK_zhemm_Astat(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k),  /* ldak * X */
                                           B(k, n),  /* ldb  * Y */
                                    zone,  C(m, n)); /* ldc  * Y */
                        }
                        else {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamConjTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m),  /* ldak * X */
                                       B(k, n),  /* ldb  * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                    }
                }
                /*
                 *  ChamLeft / ChamUpper
                 */
                else {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = C->get_blkdim( C, k, DIM_m, C->m );

                        if (k < m) {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamConjTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m),  /* ldak * X */
                                       B(k, n),  /* ldb  * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                        else if (k == m) {
                            INSERT_TASK_zhemm_Astat(
                                options,
                                side, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k),  /* ldak * K */
                                       B(k, n),  /* ldb  * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                        else {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k),  /* lda * K */
                                       B(k, n),  /* ldb * Y */
                                zone,  C(m, n)); /* ldc * Y */
                        }
                    }
                }
            }
            /*
             *  ChamRight / ChamLower
             */
            else {
                if (uplo == ChamLower) {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = C->get_blkdim( C, k, DIM_n, C->n );

                        if (k < n) {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamConjTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb * K */
                                       A(n, k),  /* lda * K */
                                zone,  C(m, n)); /* ldc * Y */
                        }
                        else if (k == n) {
                            INSERT_TASK_zhemm_Astat(
                                options,
                                side, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k),  /* ldak * Y */
                                       B(m, k),  /* ldb  * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                        else {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb  * K */
                                       A(k, n),  /* ldak * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                    }
                }
                /*
                 *  ChamRight / ChamUpper
                 */
                else {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = C->get_blkdim( C, k, DIM_n, C->n );

                        if (k < n) {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb  * K */
                                       A(k, n),  /* ldak * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                        else if (k == n) {
                            INSERT_TASK_zhemm_Astat(
                                options,
                                side, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k),  /* ldak * Y */
                                       B(m, k),  /* ldb  * Y */
                                zone,  C(m, n)); /* ldc  * Y */
                        }
                        else {
                            INSERT_TASK_zgemm_Astat(
                                options,
                                ChamNoTrans, ChamConjTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb * K */
                                       A(n, k),  /* lda * K */
                                zone,  C(m, n)); /* ldc * Y */
                        }
                    }
                }
            }

            RUNTIME_zgersum_submit_tree( options, C(m, n) );
            chameleon_data_flush( sequence, C(m, n), request->flush );
        }
    }
    options->forcesub = 0;
    free( reduceC );

    (void)chamctxt;
}


/**
 *  Parallel tile hermitian matrix-matrix multiplication.
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzhemm_summa_left( CHAM_context_t *chamctxt, cham_uplo_t uplo,
                             CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                             CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                             CHAM_desc_t *WA, CHAM_desc_t *WB,
                             RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;
    RUNTIME_request_t  *request  = options->request;
    cham_trans_t transA;
    int m, n, k, KT, lp, lq;
    int tempmm, tempnn, tempkk;
    int lookahead, myp, myq;
    int P, Q;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    lookahead = chamctxt->lookahead;
    KT  = A->nt;
    P   = chameleon_desc_datadist_get_iparam(C, 0);
    Q   = chameleon_desc_datadist_get_iparam(C, 1);
    myp = C->myrank / Q;
    myq = C->myrank % Q;

    for (k = 0; k < KT; k++ ) {
        lp = (k % lookahead) * P;
        lq = (k % lookahead) * Q;
        tempkk = A->get_blkdim( A, k, DIM_n, A->n );
        zbeta = k == 0 ? beta : zone;

        /* Transfert ownership of the k column of A or B */
        for (m = 0; m < C->mt; m ++ ) {
            int Am, Ak;

            if ( (( uplo == ChamUpper ) && ( m > k )) ||
                 (( uplo == ChamLower ) && ( m < k )) )
            {
                    /* Let's take A( k, m ) */
                Am = k;
                Ak = m;
            }
            else {
                /* Let's take A( m, k ) */
                Am = m;
                Ak = k;
            }

            chameleon_pzbcast_tile( ChamRowwise, ChamBcastRing,
                                    A( Am, Ak ), WA( m, lq ), options );
            chameleon_data_flush( sequence, A( Am, Ak ), request->flush );
        }

        /* Transfert ownership of the k row of B, or A */
        for (n = 0; n < C->nt; n++) {
            chameleon_pzbcast_tile( ChamColumnwise, ChamBcastRing,
                                    B( k, n ), WB( lp, n ), options );
            chameleon_data_flush( sequence, B( k, n ), request->flush );
        }

        /* Perform the update of this iteration */
        for (m = myp; m < C->mt; m+=P) {
            tempmm = C->get_blkdim( C, m, DIM_m, C->m );

            if ( k == m ) {
                for (n = myq; n < C->nt; n+=Q) {
                    tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                    INSERT_TASK_zhemm(
                        options, ChamLeft, uplo,
                        tempmm, tempnn, A->mb,
                        alpha, WA( m,        myq + lq ),
                               WB( myp + lp, n        ),
                        zbeta, C(  m,        n        ) );
                }
            }
            else {
                if ( (( uplo == ChamUpper ) && ( m > k )) ||
                     (( uplo == ChamLower ) && ( m < k )) )
                {
                    transA = ChamConjTrans;
                }
                else {
                    transA = ChamNoTrans;
                }

                for (n = myq; n < C->nt; n+=Q) {
                    tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                    INSERT_TASK_zgemm(
                        options, transA, ChamNoTrans,
                        tempmm, tempnn, tempkk, A->mb,
                        alpha, WA( m,        myq + lq ),
                               WB( myp + lp, n        ),
                        zbeta, C(  m,        n        ) );
                }
            }
        }
    }
}

/**
 *  Parallel tile hermitian matrix-matrix multiplication.
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzhemm_summa_right( CHAM_context_t *chamctxt, cham_uplo_t uplo,
                              CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                              CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                              CHAM_desc_t *WA, CHAM_desc_t *WB,
                              RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;
    RUNTIME_request_t  *request  = options->request;
    cham_trans_t transA;
    int m, n, k, KT, lp, lq;
    int tempmm, tempnn, tempkk;
    int lookahead, myp, myq;
    int P, Q;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    lookahead = chamctxt->lookahead;
    KT  = A->mt;
    P   = chameleon_desc_datadist_get_iparam(C, 0);
    Q   = chameleon_desc_datadist_get_iparam(C, 1);
    myp = C->myrank / Q;
    myq = C->myrank % Q;

    for (k = 0; k < KT; k++ ) {
        lp = (k % lookahead) * P;
        lq = (k % lookahead) * Q;
        tempkk = A->get_blkdim( A, k, DIM_m, A->m );
        zbeta = k == 0 ? beta : zone;

        /* Transfert ownership of the k column of A or B */
        for (m = 0; m < C->mt; m++ ) {

            chameleon_pzbcast_tile( ChamRowwise, ChamBcastRing,
                                    B( m, k ), WA( m, lq ), options );
            chameleon_data_flush( sequence, B( m, k ), request->flush );
        }

        /* Transfert ownership of the k row of B, or A */
        for (n = 0; n < C->nt; n++) {
            int Ak, An;

            if ( (( uplo == ChamUpper ) && ( n < k )) ||
                 (( uplo == ChamLower ) && ( n > k )) )
            {
                Ak = n;
                An = k;
            }
            else
            {
                Ak = k;
                An = n;
            }

            chameleon_pzbcast_tile( ChamColumnwise, ChamBcastRing,
                                    A( Ak, An ), WB( lp, n ), options );
            chameleon_data_flush( sequence, A( Ak, An ), request->flush );
        }

        /* Perform the update of this iteration */
        for (n = myq; n < C->nt; n+=Q) {
            tempnn = C->get_blkdim( C, n, DIM_n, C->n );

            if ( k == n ) {
                for (m = myp; m < C->mt; m+=P) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                    /* A has been stored in WA or WB for the summa ring */
                    INSERT_TASK_zhemm(
                        options, ChamRight, uplo,
                        tempmm, tempnn, A->mb,
                        alpha, WB( myp + lp, n        ),
                               WA( m,        myq + lq ),
                        zbeta, C(  m,        n        ) );
                }
            }
            else {
                if ( (( uplo == ChamUpper ) && ( n < k )) ||
                     (( uplo == ChamLower ) && ( n > k )) )
                {
                    transA = ChamConjTrans;
                }
                else {
                    transA = ChamNoTrans;
                }

                for (m = myp; m < C->mt; m+=P) {
                    tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                    INSERT_TASK_zgemm(
                        options, ChamNoTrans, transA,
                        tempmm, tempnn, tempkk, A->mb,
                        alpha, WA( m,        myq + lq ),
                               WB( myp + lp, n        ),
                        zbeta, C(  m,        n        ) );
                }
            }
        }
    }
}

/**
 *  Parallel tile hermitian matrix-matrix multiplication.
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzhemm_summa( CHAM_context_t *chamctxt, cham_side_t side, cham_uplo_t uplo,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        CHAM_desc_t *WA, CHAM_desc_t *WB,
                        RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;

    if (side == ChamLeft) {
        chameleon_pzhemm_summa_left( chamctxt, uplo, alpha, A, B, beta, C,
                                     WA, WB, options );
    }
    else {
        chameleon_pzhemm_summa_right( chamctxt, uplo, alpha, A, B, beta, C,
                                      WA, WB, options );
    }

    CHAMELEON_Desc_Flush( WA, sequence );
    CHAMELEON_Desc_Flush( WB, sequence );
}

/**
 *  Parallel tile hermitian matrix-matrix multiplication.
 *  Generic algorithm for any data distribution.
 */
static inline void
chameleon_pzhemm_generic( CHAM_context_t *chamctxt, cham_side_t side, cham_uplo_t uplo,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                          RUNTIME_option_t *options )
{
    int k, m, n;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    for(m = 0; m < C->mt; m++) {
        tempmm = C->get_blkdim( C, m, DIM_m, C->m );
        for(n = 0; n < C->nt; n++) {
            tempnn = C->get_blkdim( C, n, DIM_n, C->n );
            /*
             *  ChamLeft / ChamLower
             */
            if (side == ChamLeft) {
                if (uplo == ChamLower) {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = C->get_blkdim( C, k, DIM_m, C->m );
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(m, k),  /* lda * K */
                                       B(k, n),  /* ldb * Y */
                                zbeta, C(m, n)); /* ldc * Y */
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zhemm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k),  /* ldak * X */
                                           B(k, n),  /* ldb  * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamConjTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(k, m),  /* ldak * X */
                                           B(k, n),  /* ldb  * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  ChamLeft / ChamUpper
                 */
                else {
                    for (k = 0; k < C->mt; k++) {
                        tempkm = C->get_blkdim( C, k, DIM_m, C->m );
                        zbeta = k == 0 ? beta : zone;
                        if (k < m) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamConjTrans, ChamNoTrans,
                                tempmm, tempnn, tempkm, A->mb,
                                alpha, A(k, m),  /* ldak * X */
                                       B(k, n),  /* ldb  * Y */
                                zbeta, C(m, n)); /* ldc  * Y */
                        }
                        else {
                            if (k == m) {
                                INSERT_TASK_zhemm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k),  /* ldak * K */
                                           B(k, n),  /* ldb  * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkm, A->mb,
                                    alpha, A(m, k),  /* lda * K */
                                           B(k, n),  /* ldb * Y */
                                    zbeta, C(m, n)); /* ldc * Y */
                            }
                        }
                    }
                }
            }
            /*
             *  ChamRight / ChamLower
             */
            else {
                if (uplo == ChamLower) {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = C->get_blkdim( C, k, DIM_n, C->n );
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamConjTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb * K */
                                       A(n, k),  /* lda * K */
                                zbeta, C(m, n)); /* ldc * Y */
                        }
                        else {
                            if (k == n) {
                                INSERT_TASK_zhemm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k),  /* ldak * Y */
                                           B(m, k),  /* ldb  * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamNoTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k),  /* ldb  * K */
                                           A(k, n),  /* ldak * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                        }
                    }
                }
                /*
                 *  ChamRight / ChamUpper
                 */
                else {
                    for (k = 0; k < C->nt; k++) {
                        tempkn = C->get_blkdim( C, k, DIM_n, C->n );
                        zbeta = k == 0 ? beta : zone;
                        if (k < n) {
                            INSERT_TASK_zgemm(
                                options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldb  * K */
                                       A(k, n),  /* ldak * Y */
                                zbeta, C(m, n)); /* ldc  * Y */
                        }
                        else {
                            if (k == n) {
                                INSERT_TASK_zhemm(
                                    options,
                                    side, uplo,
                                    tempmm, tempnn, A->mb,
                                    alpha, A(k, k),  /* ldak * Y */
                                           B(m, k),  /* ldb  * Y */
                                    zbeta, C(m, n)); /* ldc  * Y */
                            }
                            else {
                                INSERT_TASK_zgemm(
                                    options,
                                    ChamNoTrans, ChamConjTrans,
                                    tempmm, tempnn, tempkn, A->mb,
                                    alpha, B(m, k),  /* ldb * K */
                                           A(n, k),  /* lda * K */
                                    zbeta, C(m, n)); /* ldc * Y */
                            }
                        }
                    }
                }
            }
        }
    }
    (void)chamctxt;
}

/**
 *  Parallel tile hermitian matrix-matrix multiplication. wrapper.
 */
void
chameleon_pzhemm( struct chameleon_pzgemm_s *ws,
                  cham_side_t side, cham_uplo_t uplo,
                  CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                  CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    cham_gemm_t alg = (ws != NULL) ? ws->alg : ChamGemmAlgGeneric;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    switch( alg ) {
    case ChamGemmAlgAuto:
    case ChamGemmAlgSummaB: /* Switch back to generic since it does not exist yet. */
    case ChamGemmAlgGeneric:
        chameleon_pzhemm_generic( chamctxt, side, uplo, alpha, A, B, beta, C, &options );
        break;

    case ChamGemmAlgSummaC:
        chameleon_pzhemm_summa( chamctxt, side, uplo, alpha, A, B, beta, C,
                                &(ws->WA), &(ws->WB), &options );
        break;

    case ChamGemmAlgSummaA:
        if ( side == ChamLeft ) {
            chameleon_pzhemm_Astat( chamctxt, side, uplo, alpha, A, B, beta, C, &options );
        }
        else {
            chameleon_pzhemm_generic( chamctxt, side, uplo, alpha, A, B, beta, C, &options );
        }
        break;
    }
    /* Mark written data for synchronization */
    C->sync = 1;

    RUNTIME_options_finalize( &options, chamctxt );
}
