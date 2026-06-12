/**
 *
 * @file quark/codelet_ztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpqrt Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static void
CORE_ztpqrt_quark( Quark *quark )
{
    int M;
    int N;
    int L;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileT;
    CHAMELEON_Complex64_t *WORK;

    quark_unpack_args_8( quark, M, N, L, ib,
                          tileA, tileB, tileT, WORK );

    TCORE_zlaset( ChamUpperLower, ib, N, 0., 0., tileT );
    TCORE_ztpqrt( M, N, L, ib,
                  tileA, tileB, tileT, WORK );
}

void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TPQRT;

    int shapeB = ( L == 0 ) ? 0 : (QUARK_REGION_U | QUARK_REGION_D);

    QUARK_Insert_Task(
        opt->quark, CORE_ztpqrt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                         &M,   VALUE,
        sizeof(int),                         &N,   VALUE,
        sizeof(int),                         &L,   VALUE,
        sizeof(int),                         &ib,  VALUE,
        sizeof(void*), RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), INOUT | QUARK_REGION_U | QUARK_REGION_D,
        sizeof(void*), RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), INOUT | shapeB | LOCALITY,
        sizeof(void*), RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), OUTPUT,
        sizeof(CHAMELEON_Complex64_t)*(ib+1)*nb,  NULL, SCRATCH,
        0);
}
