/**
 *
 * @file starpu/codelet_ipiv.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to work with ipiv array
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "chameleon_starpu_internal.h"

struct cl_ipiv_to_perm_args_s {
    int withidx;
    int m0;
    int m;
    int k;
    int K1;
    int K2;
    int mt;
};

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_desc_pivot_t *pivot, int k, int h, int rank )
{
    starpu_data_handle_t prevpiv = RUNTIME_pivot_getaddr( pivot, rank, h-1 );

#if defined(HAVE_STARPU_MPI_REDUX) && defined(CHAMELEON_USE_MPI)
#if !defined(HAVE_STARPU_MPI_REDUX_WRAPUP)
    starpu_data_handle_t nextpiv = RUNTIME_pivot_getaddr( pivot, rank, h );
    if ( h < pivot->n ) {
        starpu_mpi_redux_data_prio_tree( options->sequence->comm, nextpiv,
                                         options->priority, 2 /* Binary tree */ );
    }
#endif
#endif

    /* Invalidate the previous pivot structure for correct initialization in later reuse */
    if ( h > 0 ) {
        starpu_data_invalidate_submit( prevpiv );
    }

    (void)options;
    (void)k;
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_ipiv_to_perm_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_ipiv_to_perm_args_s *clargs = (struct cl_ipiv_to_perm_args_s *)cl_arg;
    int *ipiv, *perm, *invp;

    ipiv = (int*)STARPU_VECTOR_GET_PTR(descr[0]);
    perm = (int*)STARPU_VECTOR_GET_PTR(descr[1]);
    invp = (int*)STARPU_VECTOR_GET_PTR(descr[2]);

    if ( clargs->withidx ) {
        int permtmp[ clargs->m ];
        int invptmp[ clargs->m ];
        CORE_ipiv_to_perm( clargs->m0, clargs->m, clargs->k,
                           clargs->K1, clargs->K2,
                           ipiv, permtmp, invptmp );
        CORE_perm_to_idx( clargs->m0, clargs->m, clargs->mt,
                          permtmp, invptmp, perm, invp );
    }
    else {
        CORE_ipiv_to_perm( clargs->m0, clargs->m, clargs->k,
                           clargs->K1, clargs->K2,
                           ipiv, perm, invp );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
* Codelet definition
*/
static struct starpu_codelet cl_ipiv_to_perm = {
    .where        = STARPU_CPU,
#if defined(CHAMELEON_SIMULATION)
    .cpu_funcs[0] = (starpu_cpu_func_t)1,
#else
    .cpu_funcs[0] = cl_ipiv_to_perm_cpu_func,
#endif
    .nbuffers     = 3,
    .model        = NULL,
    .name         = "ipiv_to_perm"
};

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2, int mt,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    struct cl_ipiv_to_perm_args_s *clargs = NULL;

    /* Handle cache */
    if ( ipivdesc->get_rankof( ipivdesc, ipivk, ipivk ) != ipivdesc->myrank ) {
        return;
    }

    clargs = malloc( sizeof(struct cl_ipiv_to_perm_args_s) );
    clargs->withidx = ipivdesc->withidx;
    clargs->m0 = m0;
    clargs->m  = m;
    clargs->k  = k;
    clargs->K1 = K1;
    clargs->K2 = K2;
    clargs->mt = mt;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ipiv_to_perm,

        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ipiv_to_perm_args_s),

        /* Task handles */
        STARPU_R, RUNTIME_ipiv_getaddr( ipivdesc, ipivk ),
        STARPU_W, RUNTIME_ipiv_getperm( ipivdesc, ipivk ),
        STARPU_W, RUNTIME_ipiv_getinvp( ipivdesc, ipivk ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS_NOCB,
        0 );
}

