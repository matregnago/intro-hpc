/**
 *
 * @file starpu/codelet_zplssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplssq StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq_cpu_func(void *descr[], void *cl_arg)
{
    cham_store_t storev;
    int M;
    int N;
    CHAM_tile_t *tileIN;
    CHAM_tile_t *tileOUT;

    starpu_codelet_unpack_args( cl_arg, &storev, &M, &N );
    tileIN  = cti_interface_get(descr[0]);
    tileOUT = cti_interface_get(descr[1]);

    assert( tileIN->format  & CHAMELEON_TILE_FULLRANK );
    assert( tileOUT->format & CHAMELEON_TILE_FULLRANK );

    CORE_zplssq( storev, M, N, tileIN->mat, tileOUT->mat );

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq, cl_zplssq_cpu_func)

void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *IN,  int INm,  int INn,
                         const CHAM_desc_t *OUT, int OUTm, int OUTn )
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(  IN,  INm,  INn  );
    CHAMELEON_ACCESS_RW( OUT, OUTm, OUTn );
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zplssq,
        STARPU_VALUE,    &storev,            sizeof(int),
        STARPU_VALUE,    &M,                 sizeof(int),
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_R,  RTBLKADDR( IN,  ChamRealDouble, INm,  INn  ),
        STARPU_RW, RTBLKADDR( OUT, ChamRealDouble, OUTm, OUTn ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zplssq ),
        0 );
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplssq2_cpu_func(void *descr[], void *cl_arg)
{
    int N;
    CHAM_tile_t *tileRESULT;

    starpu_codelet_unpack_args(cl_arg, &N);
    tileRESULT = cti_interface_get(descr[0]);

    assert( tileRESULT->format  & CHAMELEON_TILE_FULLRANK );

    CORE_zplssq2( N, tileRESULT->mat );

    (void)cl_arg;
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplssq2, cl_zplssq2_cpu_func)

void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn )
{

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( RESULT, RESULTm, RESULTn );
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_zplssq2,
        STARPU_VALUE,    &N,                 sizeof(int),
        STARPU_RW, RTBLKADDR(RESULT, ChamRealDouble, RESULTm, RESULTn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zplssq2 ),
        0 );
}
