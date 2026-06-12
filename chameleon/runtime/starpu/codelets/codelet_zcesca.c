/**
 *
 * @file starpu/codelet_zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca StarPU codelet
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zcesca_args_s {
    int center;
    int scale;
    cham_store_t axis;
    int m;
    int n;
    int mt;
    int nt;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zcesca_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zcesca_args_s *clargs = (struct cl_zcesca_args_s *)cl_arg;
    CHAM_tile_t *Gi;
    CHAM_tile_t *Gj;
    CHAM_tile_t *G;
    CHAM_tile_t *Di;
    CHAM_tile_t *Dj;
    CHAM_tile_t *A;

    Gi = cti_interface_get(descr[0]);
    Gj = cti_interface_get(descr[1]);
    G  = cti_interface_get(descr[2]);
    Di = cti_interface_get(descr[3]);
    Dj = cti_interface_get(descr[4]);
    A  = cti_interface_get(descr[5]);

    TCORE_zcesca( clargs->center, clargs->scale, clargs->axis,
                  clargs->m, clargs->n, clargs->mt, clargs->nt,
                  Gi, Gj, G, Di, Dj, A );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zcesca, cl_zcesca_cpu_func)

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
    struct cl_zcesca_args_s *clargs = NULL;
    int exec = 0;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(Gi, Gim, Gin);
    CHAMELEON_ACCESS_R(Gj, Gjm, Gjn);
    CHAMELEON_ACCESS_R(G, Gm, Gn);
    CHAMELEON_ACCESS_R(Di, Dim, Din);
    CHAMELEON_ACCESS_R(Dj, Djm, Djn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zcesca_args_s ) );
        clargs->center = center;
        clargs->scale  = scale;
        clargs->axis   = axis;
        clargs->m      = m;
        clargs->n      = n;
        clargs->mt     = mt;
        clargs->nt     = nt;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zcesca,

        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zcesca_args_s),

        /* Task handles */
        STARPU_R,  RTBLKADDR(Gi, ChamComplexDouble, Gim, Gin),
        STARPU_R,  RTBLKADDR(Gj, ChamComplexDouble, Gjm, Gjn),
        STARPU_R,  RTBLKADDR(G,  ChamComplexDouble, Gm,  Gn ),
        STARPU_R,  RTBLKADDR(Di, ChamRealDouble,    Dim, Din),
        STARPU_R,  RTBLKADDR(Dj, ChamRealDouble,    Djm, Djn),
        STARPU_RW, RTBLKADDR(A,  ChamComplexDouble, Am,  An ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zcesca ),
        0 );
}
