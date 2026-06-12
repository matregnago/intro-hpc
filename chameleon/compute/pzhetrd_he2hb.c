/**
 *
 * @file pzhetrd_he2hb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhetrd_he2hb parallel algorithm
 *
 * @version 1.4.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m, n) A,  m,  n
#define T(m, n) T,  m,  n
#define D(k)  &D, (k)-1, 0

#define AT(k) &AT, k, 0

#if defined(CHAMELEON_COPY_DIAG)
#define E(m, n) E, m, 0
#else
#define E(m, n) A, m, n
#endif

/**
 *  Parallel tile BAND Tridiagonal Reduction - dynamic scheduler
 */
void chameleon_pzhetrd_he2hb( cham_uplo_t uplo,
                              CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *E,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t D;
    CHAM_desc_t AT;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, j;
    int tempkm, tempkn, tempmm, tempnn, tempjj;
    int ib;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }

    RUNTIME_options_init(&options, chamctxt, sequence, request);
    ib = CHAMELEON_IB;

    /*
     * zgeqrt        = A->nb * (ib+1)
     * zunmqr        = A->nb * ib
     * ztsqrt        = A->nb * (ib+1)
     * ztsmqr        = A->nb * ib
     * zherfb        = A->nb * ib
     * ztsmqr_hetra1 = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr =     A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     * zherfb =     A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Copy of the diagonal tiles to keep the general version of the tile all along the computation */
    chameleon_ddesc_alloc_diag( chamctxt, &D, "HETRD_HE2HB_D", A->mb, A->m, A->n,
                                chameleon_desc_datadist_get_iparam(A, 0),
                                chameleon_desc_datadist_get_iparam(A, 1) );

    chameleon_desc_init_2dlap( chamctxt, &AT, "HETRD_HE2HB_AT", ChamComplexDouble,
                               A->mb, A->nb,
                               chameleon_min(A->mt, A->nt) * A->mb, A->nb, 1, 1 );

    /* Let's extract the diagonal in a temporary copy that contains A and A' */
    for (k = 1; k < A->nt; k++){
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );

        INSERT_TASK_zhe2ge( &options,
                            uplo, tempkn, tempkn, A->mb,
                            A(k, k), D(k) );
    }

    if (uplo == ChamLower) {
       for (k = 0; k < A->nt-1; k++){
           RUNTIME_iteration_push(chamctxt, k);

           tempkm = A->get_blkdim( A, k+1, DIM_m, A->m );
           tempkn = A->get_blkdim( A, k,   DIM_n, A->n );

           INSERT_TASK_zgeqrt(
               &options,
               tempkm, tempkn, ib, A->nb,
               A(k+1, k),
               T(k+1, k));

#if defined(CHAMELEON_COPY_DIAG)
           INSERT_TASK_zlacpy(
               &options,
               ChamLower, tempkm, tempkn,
               A(k+1, k),
               E(k+1, k) );
#if defined(CHAMELEON_USE_CUDA)
           INSERT_TASK_zlaset(
               &options,
               ChamUpper, tempkm, tempkn,
               0., 1.,
               E(k+1, k) );
#endif
#endif

           /* LEFT and RIGHT on the symmetric diagonal block */
           INSERT_TASK_zherfb(
               &options,
               ChamLower,
               tempkm, tempkm, ib, A->nb,
               E(k+1, k),
               T(k+1, k),
               D(k+1));

           /* RIGHT on the remaining tiles until the bottom */
           for (m = k+2; m < A->mt ; m++) {
               tempmm = A->get_blkdim( A, m, DIM_m, A->m );
               INSERT_TASK_zunmqr(
                   &options,
                   ChamRight, ChamNoTrans,
                   tempmm, A->nb, tempkm, ib, A->nb,
                   E(k+1, k),
                   T(k+1, k),
                   A(m,   k+1));
           }

           for (m = k+2; m < A->mt; m++) {
               tempmm = A->get_blkdim( A, m, DIM_m, A->m );

               options.priority = 1;
               INSERT_TASK_ztsqrt(
                   &options,
                   tempmm, A->nb, ib, A->nb,
                   A(k+1, k),
                   A(m  , k),
                   T(m  , k));
               options.priority = 0;

               /* LEFT */
               for (i = k+2; i < m; i++) {
                   INSERT_TASK_ztsmqr_hetra1(
                       &options,
                       ChamLeft, ChamConjTrans,
                       A->mb, A->nb, tempmm, A->nb, A->nb, ib, A->nb,
                       A(i, k+1),
                       A(m,   i),
                       A(m,   k),
                       T(m,   k));
               }

               /* RIGHT */
               for (j = m+1; j < A->mt ; j++) {
                   tempjj = A->get_blkdim( A, j, DIM_m, A->m );
                   INSERT_TASK_ztsmqr(
                       &options,
                       ChamRight, ChamNoTrans,
                       tempjj, A->nb, tempjj, tempmm, A->nb, ib, A->nb,
                       A(j, k+1),
                       A(j,   m),
                       A(m,   k),
                       T(m,   k));
               }

               /* LEFT->RIGHT */
               options.priority = 1;
               /**
                * Compute the update of the following:
                *
                *     A1   A2'
                *     A2   A3
                *
                * with A1 and A3 two diagonal tiles. This is the tsmqr_corner
                * from plasma split in 4 tasks
                */
               /*  Copy the transpose of A2 (m, k+1): AT(k) <- A2' = A2(k+1, m) */
               INSERT_TASK_zlatro(
                   &options,
                   ChamUpperLower, ChamConjTrans,
                   tempmm, A->nb, A->nb,
                   A(m, k+1),
                   AT(m));

               /*  Left application on |A1| */
               /*                      |A2| */
               INSERT_TASK_ztsmqr(
                   &options,
                   ChamLeft, ChamConjTrans,
                   A->mb, A->nb, tempmm, A->nb, A->nb, ib, A->nb,
                   D(k+1),
                   A(m, k+1),
                   A(m,   k),
                   T(m,   k));

               /*  Left application on | A2'| */
               /*                      | A3 | */
               INSERT_TASK_ztsmqr(
                   &options,
                   ChamLeft, ChamConjTrans,
                   A->mb, tempmm, tempmm, tempmm, A->nb, ib, A->nb,
                   AT(m),
                   D(m),
                   A(m,  k),
                   T(m,  k));

               /*  Right application on | A1 A2' | */
               INSERT_TASK_ztsmqr(
                   &options,
                   ChamRight, ChamNoTrans,
                   A->mb, A->nb, A->mb, tempmm, A->nb, ib, A->nb,
                   D(k+1),
                   AT(m) ,
                   A(m,   k),
                   T(m,   k));

               /*  Right application on | A2 A3 | */
               INSERT_TASK_ztsmqr(
                   &options,
                   ChamRight, ChamNoTrans,
                   tempmm, A->nb, tempmm, tempmm, A->nb, ib, A->nb,
                   A(m, k+1),
                   D(m),
                   A(m,   k),
                   T(m,   k));
               options.priority = 0;
           }

           RUNTIME_iteration_pop(chamctxt);
       }
    }
    else {
       for (k = 0; k < A->nt-1; k++){
           RUNTIME_iteration_push(chamctxt, k);

           tempkm = A->get_blkdim( A, k,   DIM_m, A->m );
           tempkn = A->get_blkdim( A, k+1, DIM_n, A->n );

           INSERT_TASK_zgelqt(
               &options,
               tempkm, tempkn, ib, A->nb,
               A(k, k+1),
               T(k, k+1));

#if defined(CHAMELEON_COPY_DIAG)
           INSERT_TASK_zlacpy(
               &options,
               ChamUpper, tempkm, tempkn,
               A(k, k+1),
               E(k, k+1) );
#if defined(CHAMELEON_USE_CUDA)
           INSERT_TASK_zlaset(
               &options,
               ChamLower, tempkm, tempkn,
               0., 1.,
               E(k, k+1) );
#endif
#endif

           /* RIGHT and LEFT on the symmetric diagonal block */
           INSERT_TASK_zherfb(
               &options,
               ChamUpper,
               tempkn, tempkn, ib, A->nb,
               E(k, k+1),
               T(k, k+1),
               D(k+1));

           /* LEFT on the remaining tiles until the left side */
           for (n = k+2; n < A->nt ; n++) {
               tempnn = A->get_blkdim( A, n, DIM_n, A->n );
               INSERT_TASK_zunmlq(
                   &options,
                   ChamLeft, ChamNoTrans,
                   A->mb, tempnn, tempkn, ib, A->nb,
                   E(k,   k+1),
                   T(k,   k+1),
                   A(k+1, n  ));
           }

           for (n = k+2; n < A->nt; n++) {
               tempnn = A->get_blkdim( A, n, DIM_n, A->n );
               options.priority = 1;
               INSERT_TASK_ztslqt(
                   &options,
                   A->mb, tempnn, ib, A->nb,
                   A(k, k+1),
                   A(k, n  ),
                   T(k, n  ));
               options.priority = 0;

               /* RIGHT */
               for (i = k+2; i < n; i++) {
                   INSERT_TASK_ztsmlq_hetra1(
                       &options,
                       ChamRight, ChamConjTrans,
                       A->mb, A->nb, A->nb, tempnn, A->nb, ib, A->nb,
                       A(k+1, i),
                       A(i,   n),
                       A(k,   n),
                       T(k,   n));
               }

               /* LEFT */
               for (j = n+1; j < A->nt ; j++) {
                   tempjj = A->get_blkdim( A, j, DIM_n, A->n );
                   INSERT_TASK_ztsmlq(
                       &options,
                       ChamLeft, ChamNoTrans,
                       A->nb, tempjj, tempnn, tempjj, A->nb, ib, A->nb,
                       A(k+1, j),
                       A(n,   j),
                       A(k,   n),
                       T(k,   n));
               }

               /* RIGHT->LEFT */
               options.priority = 1;
               /**
                * Compute the update of the following:
                *
                *     A1   A2
                *     A2'  A3
                *
                * with A1 and A3 two diagonal tiles. This is the tsmqr_corner
                * from plasma split in 4 tasks
                */
               /*  Copy the transpose of A2: AT(k) <- A2' */
               INSERT_TASK_zlatro(
                   &options,
                   ChamUpperLower, ChamConjTrans,
                   A->mb, tempnn, A->nb,
                   A(k+1, n),
                   AT(n) );

               /*  Right application on | A1 A2 | */
               INSERT_TASK_ztsmlq(
                   &options,
                   ChamRight, ChamConjTrans,
                   A->mb, A->nb, A->mb, tempnn, A->nb, ib, A->nb,
                   D(k+1),
                   A(k+1, n),
                   A(k,   n),
                   T(k,   n));

               /*  Right application on | A2' A3 | */
               INSERT_TASK_ztsmlq(
                   &options,
                   ChamRight, ChamConjTrans,
                   tempnn, A->nb, tempnn, tempnn, A->nb, ib, A->nb,
                   AT(n),
                   D(n),
                   A(k,  n),
                   T(k,  n));

               /*  Left application on |A1 | */
               /*                      |A2'| */
               INSERT_TASK_ztsmlq(
                   &options,
                   ChamLeft, ChamNoTrans,
                   A->mb, A->nb, tempnn, A->nb, A->nb, ib, A->nb,
                   D(k+1),
                   AT(n),
                   A(k, n),
                   T(k, n));

               /*  Left application on | A2 | */
               /*                      | A3 | */
               INSERT_TASK_ztsmlq(
                   &options,
                   ChamLeft, ChamNoTrans,
                   A->mb, tempnn, tempnn, tempnn, A->nb, ib, A->nb,
                   A(k+1, n),
                   D(n),
                   A(k,   n),
                   T(k,   n));
           }
           options.priority = 0;

           RUNTIME_iteration_pop(chamctxt);
       }
    }

    /* Copy-back into A */
    for (k = 1; k < A->nt; k++){
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        INSERT_TASK_zlacpy( &options,
                            uplo, tempkn, tempkn,
                            D(k), A(k, k));
    }


    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);

    CHAMELEON_Sequence_Wait(sequence);
    chameleon_desc_destroy( &D );
    chameleon_desc_destroy( &AT );

    (void)E;
}
