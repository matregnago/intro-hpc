/**
 *
 * @file quark/codelet_zplgsy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplgsy Quark codelet
 *
 * @version 1.4.0
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
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zplgsy_quark(Quark *quark)
{
    CHAMELEON_Complex64_t bump;
    int m;
    int n;
    CHAM_tile_t *tileA;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_8( quark, bump, m, n, tileA, bigM, m0, n0, seed );
    TCORE_zplgsy( bump, m, n, tileA, bigM, m0, n0, seed );
}

void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                        CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLGSY;
    QUARK_Insert_Task(opt->quark, CORE_zplgsy_quark, (Quark_Task_Flags*)opt,
        sizeof(CHAMELEON_Complex64_t),       &bump, VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),         OUTPUT,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}
