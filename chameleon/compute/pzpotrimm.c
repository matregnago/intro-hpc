/**
 *
 * @file pzpotrimm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrimm parallel algorithm
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Ali M Charara
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-12-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void chameleon_pzpotrimm(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempmm, tempnn, tempkn;

    CHAMELEON_Complex64_t alpha = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t beta  = (CHAMELEON_Complex64_t) 0.0;
    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;


    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );

            INSERT_TASK_zpotrf(
                &options,
                ChamLower, tempkm, A->mb,
                A(k, k), A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, ChamLower, ChamConjTrans, ChamNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k),
                          A(m, k));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_zherk(
                    &options,
                    ChamLower, ChamNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k),
                     1.0, A(n, n));

                for (m = n+1; m < A->mt; m++) {
                    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k),
                               A(n, k),
                        zone,  A(m, n));
                }
                chameleon_data_flush( sequence, A(n, k), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, A->nt + k);

            tempkn = A->get_blkdim( A, k, DIM_n, A->n );
            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, ChamNonUnit,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k),
                           A(m, k));
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                for (n = 0; n < k; n++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k),
                              A(k, n),
                        zone, A(m, n));
                }
                chameleon_data_flush( sequence, A(m, k), request->flush );
            }
            for (n = 0; n < k; n++) {
                chameleon_data_flush( sequence, A(k, n), request->flush );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, ChamNonUnit,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k),
                          A(k, n));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_ztrtri(
                &options,
                uplo, ChamNonUnit,
                tempkn, A->mb,
                A(k, k), A->nb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, 2*A->nt + k);

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            for(n = 0; n < k; n++) {
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n),
                    1.0, A(n, n));

                for(m = n+1; m < k; m++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m),
                             A(k, n),
                        1.0, A(m, n));
                }
            }
            for (n = 0; n < k; n++) {
                chameleon_data_flush( sequence, A(k, n), request->flush );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamLeft, uplo, ChamConjTrans, ChamNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k),
                         A(k, n));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k));

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZSYMM Right / Lower
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(chamctxt, 3*A->nt + k);

            tempkn = C->get_blkdim( C, k, DIM_n, C->n );
            zbeta = k == 0 ? beta : zone;

            for (m = 0; m < C->mt; m++) {
                tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                for (n = 0; n < C->nt; n++) {
                    tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                    if (k < n) {
                       INSERT_TASK_zgemm(
                           &options,
                           ChamNoTrans, ChamTrans,
                           tempmm, tempnn, tempkn, A->mb,
                           alpha, B(m, k),  /* ldbm * K */
                                  A(n, k),  /* ldan * K */
                           zbeta, C(m, n)); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                           INSERT_TASK_zsymm(
                               &options,
                               ChamRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k),  /* ldak * Y */
                                      B(m, k),  /* ldbm * Y */
                               zbeta, C(m, n)); /* ldcm * Y */
                        }
                        else {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldbm * K */
                                       A(k, n),  /* ldak * Y */
                                zbeta, C(m, n)); /* ldcm * Y */
                        }
                    }
                }
                chameleon_data_flush( sequence, B(m, k), request->flush );
            }
            for (n = 0; n <= k; n++) {
                chameleon_data_flush( sequence, A(k, n), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /*
     *  ChamUpper
     */
    else {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = A->get_blkdim( A, k, DIM_n, A->n );
            INSERT_TASK_zpotrf(
                &options,
                ChamUpper,
                tempkm, A->mb,
                A(k, k), A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, ChamUpper, ChamConjTrans, ChamNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k),
                          A(k, n));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );

                INSERT_TASK_zherk(
                    &options,
                    ChamUpper, ChamConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m),
                     1.0, A(m, m));

                for (n = m+1; n < A->nt; n++) {
                    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m),
                               A(k, n),
                        zone,  A(m, n));
                }
                chameleon_data_flush( sequence, A(k, m), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, A->nt + k);

            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, uplo, ChamNoTrans, ChamNonUnit,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k),
                           A(k, n));
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                for (m = 0; m < k; m++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k),
                              A(k, n),
                        zone, A(m, n));
                }
                chameleon_data_flush( sequence, A(k, n), request->flush );
            }
            for (m = 0; m < k; m++) {
                chameleon_data_flush( sequence, A(m, k), request->flush );
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, uplo, ChamNoTrans, ChamNonUnit,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k),
                          A(m, k));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_ztrtri(
                &options,
                uplo, ChamNonUnit,
                tempkm, A->mb,
                A(k, k), A->mb*k);

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, 2*A->nt + k);

            tempkn = A->get_blkdim( A, k, DIM_n, A->n );

            for (m = 0; m < k; m++) {
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k),
                    1.0, A(m, m));

                for (n = m+1; n < k; n++){
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k),
                             A(n, k),
                        1.0, A(m, n));
                }
            }
            for (m = 0; m < k; m++) {
                chameleon_data_flush( sequence, A(m, k), request->flush );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamRight, uplo, ChamConjTrans, ChamNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k),
                         A(m, k));
            }
            chameleon_data_flush( sequence, A(k, k), request->flush );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k));

            RUNTIME_iteration_pop(chamctxt);
        }
        /*
         *  ZSYMM Right / Upper
         */
        for (k = 0; k < C->nt; k++) {
            RUNTIME_iteration_push(chamctxt, 3*A->nt + k);

            tempkn = C->get_blkdim( C, k, DIM_n, C->n );
            zbeta = k == 0 ? beta : zone;

            for (m = 0; m < C->mt; m++) {
                tempmm = C->get_blkdim( C, m, DIM_m, C->m );

                for (n = 0; n < C->nt; n++) {
                    tempnn = C->get_blkdim( C, n, DIM_n, C->n );

                    if (k < n) {
                        INSERT_TASK_zgemm(
                            &options,
                            ChamNoTrans, ChamNoTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, B(m, k),  /* ldbm * K */
                                   A(k, n),  /* ldak * Y */
                            zbeta, C(m, n)); /* ldcm * Y */
                    }
                    else {
                        if (k == n) {
                            INSERT_TASK_zsymm(
                                &options,
                                ChamRight, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k),  /* ldak * Y */
                                       B(m, k),  /* ldbm * Y */
                                zbeta, C(m, n)); /* ldcm * Y */
                        }
                        else {
                            INSERT_TASK_zgemm(
                                &options,
                                ChamNoTrans, ChamTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k),  /* ldbm * K */
                                       A(n, k),  /* ldan * K */
                                zbeta, C(m, n)); /* ldcm * Y */
                        }
                    }
                }
                chameleon_data_flush( sequence, B(m, k), request->flush );
            }
            for (m = 0; m <= k; m++) {
                chameleon_data_flush( sequence, A(m, k), request->flush );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /* Mark written data for synchronization */
    A->sync = 1;
    C->sync = 1;

    RUNTIME_options_finalize(&options, chamctxt);
}
