/**
 *
 * @file quark/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zplssq_quark(Quark *quark)
{
    cham_store_t storev;
    int M;
    int N;
    CHAM_tile_t *tileIN;
    CHAM_tile_t *tileOUT;

    quark_unpack_args_5( quark, storev, M, N, tileIN, tileOUT );

    assert( tileIN->format  & CHAMELEON_TILE_FULLRANK );
    assert( tileOUT->format & CHAMELEON_TILE_FULLRANK );
    CORE_zplssq( storev, M, N, tileIN->mat, tileOUT->mat );
}

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *IN,  int INm,  int INn,
                         const CHAM_desc_t *OUT, int OUTm, int OUTn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLSSQ;
    QUARK_Insert_Task(opt->quark, CORE_zplssq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),            &storev,    VALUE,
        sizeof(int),            &M,         VALUE,
        sizeof(int),            &N,         VALUE,
        sizeof(void*), RTBLKADDR(IN,  double, INm,  INn),  INPUT,
        sizeof(void*), RTBLKADDR(OUT, double, OUTm, OUTn), INOUT,
        0);
}


void CORE_zplssq2_quark(Quark *quark)
{
    int N;
    CHAM_tile_t *tileRESULT;

    quark_unpack_args_2( quark, N, tileRESULT );

    assert( tileRESULT->format & CHAMELEON_TILE_FULLRANK );
    CORE_zplssq2(N, tileRESULT->mat);
}

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLSSQ2;
    QUARK_Insert_Task(opt->quark, CORE_zplssq2_quark, (Quark_Task_Flags*)opt,
        sizeof(int),        &N,         VALUE,
        sizeof(void*), RTBLKADDR(RESULT, double, RESULTm, RESULTn), INOUT,
        0);
}
