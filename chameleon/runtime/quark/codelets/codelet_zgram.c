/**
 *
 * @file quark/codelet_zgram.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram Quark codelet
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

void CORE_zgram_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int m, n, mt, nt;
    CHAM_tile_t *Di;
    CHAM_tile_t *Dj;
    CHAM_tile_t *D;
    CHAM_tile_t *A;

    quark_unpack_args_9(quark, uplo, m, n, mt, nt, Di, Dj, D, A );
    TCORE_zgram( uplo, m, n, mt, nt, Di, Dj, D, A );
}

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GRAM;
    QUARK_Insert_Task(opt->quark, CORE_zgram_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),             &uplo,      VALUE,
                      sizeof(int),             &m,         VALUE,
                      sizeof(int),             &n,         VALUE,
                      sizeof(int),             &mt,        VALUE,
                      sizeof(int),             &nt,        VALUE,
                      sizeof(void*), RTBLKADDR(Di, double, Dim, Din), INPUT,
                      sizeof(void*), RTBLKADDR(Dj, double, Djm, Djn), INPUT,
                      sizeof(void*), RTBLKADDR(D,  double, Dm,  Dn ), INPUT,
                      sizeof(void*), RTBLKADDR(A,  double, Am,  An ), INOUT,
                      0);
}
