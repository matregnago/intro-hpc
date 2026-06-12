/**
 *
 * @file starpu/runtime_context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU context routines
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Terry Cojean
 * @author Loris Lucido
 * @date 2024-09-17
 *
 */
#include "chameleon_starpu_internal.h"

/**
 * @brief Store the status of some flags to knwo when enable/disable them
 */
static int context_starpu_flags = 0;

#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 3))
/* Defined by StarPU as external function */
#else
#if ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
int _starpu_is_initialized(void);
#define starpu_is_initialized() _starpu_is_initialized()
#else
#define starpu_is_initialized() 0
#endif
#endif

/**
 *  Create new context
 */
void RUNTIME_context_create( CHAM_context_t *chamctxt )
{
    chamctxt->scheduler = RUNTIME_SCHED_STARPU;

    if ( !starpu_is_initialized() ) {
        CHAM_context_starpu_t *sched_opt = malloc( sizeof(CHAM_context_starpu_t) );

        sched_opt->pw_config = NULL;
        starpu_conf_init( &(sched_opt->starpu_conf) );

        chamctxt->schedopt = sched_opt;
    }
    else {
        chamctxt->schedopt = NULL;
    }

    return;
}

/**
 *  Clean the context
 */
void RUNTIME_context_destroy( CHAM_context_t *chamctxt )
{
    /* StarPU was already initialized by an external library */
    if ( chamctxt->schedopt ) {
        free( chamctxt->schedopt );
    }
    return;
}

/**
 *
 */
void RUNTIME_enable( void *runtime_ctxt, int lever )
{
    switch (lever)
    {
    case CHAMELEON_DAG:
        fprintf(stderr, "StarPU is providing DAG generation through tracing support (CHAMELEON_GENERATE_TRACE)\n");
        break;
    case CHAMELEON_GENERATE_STATS:
        context_starpu_flags |= (1 << CHAMELEON_GENERATE_STATS);
        starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
        break;
    case CHAMELEON_GENERATE_TRACE:
        context_starpu_flags |= (1 << CHAMELEON_GENERATE_TRACE);
        starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
        break;
    case CHAMELEON_BOUND:
        starpu_bound_start(0, 0);
        break;
    default:
        return;
    }

    (void)runtime_ctxt;
    return;
}

/**
 *
 */
void RUNTIME_disable( void *runtime_ctxt, int lever )
{
    switch (lever)
    {
    case CHAMELEON_DAG:
        fprintf(stderr, "StarPU is providing DAG generation through tracing support (CHAMELEON_GENERATE_TRACE)\n");
        break;
    case CHAMELEON_GENERATE_TRACE:
        context_starpu_flags |= ~(1 << CHAMELEON_GENERATE_TRACE);
        if ( !context_starpu_flags ) {
            starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
        }
        break;
    case CHAMELEON_GENERATE_STATS:
        context_starpu_flags |= ~(1 << CHAMELEON_GENERATE_STATS);
        if ( !context_starpu_flags ) {
            starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
        }
        break;
    case CHAMELEON_BOUND:
        starpu_bound_stop();
        break;
    default:
        return;
    }

    (void)runtime_ctxt;
    return;
}
