/**
 *
 * @file core_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zhe2ge CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 */
void CORE_zhe2ge(cham_uplo_t uplo, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 CHAMELEON_Complex64_t *B, int LDB)
{
    const CHAMELEON_Complex64_t *Aptr;
    CHAMELEON_Complex64_t *Bptr, *BTptr;
    int i, j;

    Aptr = A;
    Bptr = B;

    if (uplo == ChamLower){
        for (j = 0; j < N; j++){
            /* Diagonal element */
            *Bptr = *Aptr;
            Bptr++; Aptr++;

            /* Outside the diagonal */
            BTptr = B + j + (j+1) * LDB;
            for (i = j+1; i < M; i++, Bptr++, Aptr++, BTptr += LDB) {
                *Bptr  = *Aptr;
                *BTptr = conj( *Aptr );
            }
            Aptr += (LDA - i + j + 1);
            Bptr += (LDB - i + j + 1);
        }
    }
    else{
        for (j = 0; j < N; j++){
            BTptr = B + j;
            for (i = 0; i < j; i++, Bptr++, Aptr++, BTptr += LDB) {
                *Bptr  = *Aptr;
                *BTptr = conj( *Aptr );
            }
            *Bptr = *A;

            Aptr += (LDA - i);
            Bptr += (LDB - i);
        }
    }
}

