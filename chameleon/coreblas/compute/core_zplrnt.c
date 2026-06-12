/**
 *
 * @file core_zplrnt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zplrnt CPU kernel
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/random.h"

#if defined(PRECISION_z) || defined(PRECISION_c)
#define NBELEM   2
#else
#define NBELEM   1
#endif

void CORE_zplrnt( int m, int n, CHAMELEON_Complex64_t *A, int lda,
                  int bigM, int m0, int n0, unsigned long long int seed )
{
    CHAMELEON_Complex64_t *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)bigM;

    for (j=0; j<n; ++j ) {
        ran = CORE_rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            *tmp = CORE_zlaran( &ran );
            tmp++;
        }
        tmp  += lda-i;
        jump += bigM;
    }
}
