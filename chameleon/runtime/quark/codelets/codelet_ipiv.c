/**
 *
 * @file quark/codelet_ipiv.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark codelets to convert pivot to permutations
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"
#include "coreblas.h"

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_desc_pivot_t *pivot, int k, int h, int rank )
{
    assert( 0 );
    (void)options;
    (void)pivot;
    (void)k;
    (void)h;
    (void)rank;
}

static inline void
CORE_ipiv_to_perm_quark( Quark *quark )
{
    int m0, m, k, K1, K2, mt;
    int *ipiv, *perm, *invp;

    quark_unpack_args_9( quark, m0, m, k, K1, K2, mt, ipiv, perm, invp );

    CORE_ipiv_to_perm( m0, m, k, K1, K2, ipiv, perm, invp );

    (void)mt;
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2, int mt,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(
        opt->quark, CORE_ipiv_to_perm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),  &m0,  VALUE,
        sizeof(int),  &m,   VALUE,
        sizeof(int),  &k,   VALUE,
        sizeof(int),  &K1,  VALUE,
        sizeof(int),  &K2,  VALUE,
        sizeof(int),  &mt,  VALUE,
        sizeof(int*), RUNTIME_ipiv_getaddr( ipivdesc, ipivk ), INPUT,
        sizeof(int*), RUNTIME_ipiv_getperm( ipivdesc, ipivk ), OUTPUT,
        sizeof(int*), RUNTIME_ipiv_getinvp( ipivdesc, ipivk ), OUTPUT,
        0 );
}
