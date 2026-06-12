/**
 *
 * @file quark/codelet_zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca Quark codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zcesca_quark(Quark *quark)
{
    int center, scale;
    cham_store_t axis;
    int m, n, mt, nt;
    CHAM_tile_t *Gi;
    CHAM_tile_t *Gj;
    CHAM_tile_t *G;
    CHAM_tile_t *Di;
    CHAM_tile_t *Dj;
    CHAM_tile_t *A;

    quark_unpack_args_13(quark, center, scale, axis, m, n, mt, nt, Gi, Gj, G, Di, Dj, A );
    TCORE_zcesca( center, scale, axis, m, n, mt, nt, Gi, Gj, G, Di, Dj, A );
}

void INSERT_TASK_zcesca( const RUNTIME_option_t *options,
                         int center, int scale, cham_store_t axis,
                         int m, int n, int mt, int nt,
                         const CHAM_desc_t *Gi, int Gim, int Gin,
                         const CHAM_desc_t *Gj, int Gjm, int Gjn,
                         const CHAM_desc_t *G, int Gm, int Gn,
                         const CHAM_desc_t *Di, int Dim, int Din,
                         const CHAM_desc_t *Dj, int Djm, int Djn,
                         CHAM_desc_t *A, int Am, int An )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_CESCA;
    QUARK_Insert_Task(opt->quark, CORE_zcesca_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),             &center,    VALUE,
                      sizeof(int),             &scale,     VALUE,
                      sizeof(int),             &axis,      VALUE,
                      sizeof(int),             &m,         VALUE,
                      sizeof(int),             &n,         VALUE,
                      sizeof(int),             &mt,        VALUE,
                      sizeof(int),             &nt,        VALUE,
                      sizeof(void*), RTBLKADDR(Gi, CHAMELEON_Complex64_t, Gim, Gin), INPUT,
                      sizeof(void*), RTBLKADDR(Gj, CHAMELEON_Complex64_t, Gjm, Gjn), INPUT,
                      sizeof(void*), RTBLKADDR(G,  CHAMELEON_Complex64_t, Gm,  Gn ), INPUT,
                      sizeof(void*), RTBLKADDR(Di, double, Dim, Din), INPUT,
                      sizeof(void*), RTBLKADDR(Dj, double, Djm, Djn), INPUT,
                      sizeof(void*), RTBLKADDR(A,  CHAMELEON_Complex64_t, Am,  An ), INOUT,
                      0);
}
