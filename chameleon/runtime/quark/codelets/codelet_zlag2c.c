/**
 *
 * @file quark/codelet_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_zc.h"
#include "coreblas/coreblas_zctile.h"

void CORE_zlag2c_quark( Quark *quark )
{
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    quark_unpack_args_6( quark, m, n, tileA, tileB, sequence, request );
    TCORE_zlag2c( m, n, tileA, tileB, &info );
    if ( (sequence->status != CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, info );
    }
}

void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAG2C;
    QUARK_Insert_Task( opt->quark, CORE_zlag2c_quark, (Quark_Task_Flags*)opt,
                       sizeof(int),                        &m,         VALUE,
                       sizeof(int),                        &n,         VALUE,
                       sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
                       sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex32_t, Bm, Bn),                 OUTPUT,
                       sizeof(RUNTIME_sequence_t*),           &(options->sequence),  VALUE,
                       sizeof(RUNTIME_request_t*),            &(options->request),   VALUE,
                       0 );

    (void)nb;
}

void CORE_clag2z_quark( Quark *quark )
{
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    quark_unpack_args_4( quark, m, n, tileA, tileB );
    TCORE_clag2z( m, n, tileA, tileB );
}

void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAG2C;

    QUARK_Insert_Task( opt->quark, CORE_clag2z_quark, (Quark_Task_Flags*)opt,
                       sizeof(int),                        &m,     VALUE,
                       sizeof(int),                        &n,     VALUE,
                       sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex32_t, Am, An),             INPUT,
                       sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             INOUT,
                       0 );
    (void)nb;
}
