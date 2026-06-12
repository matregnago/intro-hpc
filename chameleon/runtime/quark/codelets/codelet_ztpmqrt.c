/**
 *
 * @file quark/codelet_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztpmqrt Quark codelet
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
CORE_ztpmqrt_quark( Quark *quark )
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t *WORK;

    quark_unpack_args_12( quark, side, trans, M, N, K, L, ib,
                          tileV, tileT, tileA, tileB, WORK );

    TCORE_ztpmqrt( side, trans, M, N, K, L, ib,
                  tileV, tileT, tileA, tileB, WORK );
}

void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const CHAM_desc_t *V, int Vm, int Vn,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TPMQRT;

    int shapeV = ( L == 0 ) ? 0 : (QUARK_REGION_U | QUARK_REGION_D);

    QUARK_Insert_Task(
        opt->quark, CORE_ztpmqrt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &M,     VALUE,
        sizeof(int),                     &N,     VALUE,
        sizeof(int),                     &K,     VALUE,
        sizeof(int),                     &L,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(void*), RTBLKADDR( V, CHAMELEON_Complex64_t, Vm, Vn ), INPUT | shapeV,
        sizeof(void*), RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), INPUT,
        sizeof(void*), RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), INOUT,
        sizeof(void*), RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), INOUT | LOCALITY,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL, SCRATCH,
        0);
}
