/**
 *
 * @file starpu/codelet_zgram.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Lucas Barros de Assis
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgram_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int m, n, mt, nt;
    CHAM_tile_t *Di;
    CHAM_tile_t *Dj;
    CHAM_tile_t *D;
    CHAM_tile_t *A;

    Di = cti_interface_get(descr[0]);
    Dj = cti_interface_get(descr[1]);
    D  = cti_interface_get(descr[2]);
    A  = cti_interface_get(descr[3]);

    starpu_codelet_unpack_args( cl_arg, &uplo, &m, &n, &mt, &nt );
    TCORE_zgram( uplo, m, n, mt, nt, Di, Dj, D, A );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgram, cl_zgram_cpu_func)

void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An )
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(Di, Dim, Din);
    CHAMELEON_ACCESS_R(Dj, Djm, Djn);
    CHAMELEON_ACCESS_R(D, Dm, Dn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zgram,
        STARPU_VALUE,    &uplo,                      sizeof(int),
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &mt,                        sizeof(int),
        STARPU_VALUE,    &nt,                        sizeof(int),
        STARPU_R,        RTBLKADDR(Di, ChamRealDouble, Dim, Din),
        STARPU_R,        RTBLKADDR(Dj, ChamRealDouble, Djm, Djn),
        STARPU_R,        RTBLKADDR(D, ChamRealDouble, Dm, Dn),
        STARPU_RW,       RTBLKADDR(A, ChamRealDouble, Am, An),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zgram ),
        0 );
}
