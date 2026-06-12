/**
 *
 * @file tcore.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon tcore EZtrace 2.0 module
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>

/* Our own version of function none, without final coma */
#define FUNCTION_NONE_ELT  {                    \
        .function_name = "",                    \
        .callback      = NULL,                  \
        .event_id      = -1,                    \
    }

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _coreblas_tcore_initialized = 0;

struct ezt_instrumented_function *PPTRACE_SYMBOL_LIST(coreblas_tcore);

#if defined(PRECISION_z)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_z);
#endif
#if defined(PRECISION_c)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_c);
#endif
#if defined(PRECISION_d)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_d);
#endif
#if defined(PRECISION_s)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_s);
#endif

#if defined(PRECISION_zc)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_zc);
#endif
#if defined(PRECISION_ds)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_tcore_ds);
#endif

static void
init_coreblas_tcore_()
{
    INSTRUMENT_FUNCTIONS(coreblas_tcore);

    if ( eztrace_autostart_enabled() ) {
        eztrace_start();
    }

    _coreblas_tcore_initialized = 1;
}

static void
finalize_coreblas_tcore()
{
    _coreblas_tcore_initialized = 0;

    eztrace_stop();
}

static void _coreblas_tcore_init( void ) __attribute__((constructor));
static void
_coreblas_tcore_init( void )
{
    struct ezt_instrumented_function *current_list_start;
    size_t size   = 0;
    size_t sizez  = 0;
    size_t sizec  = 0;
    size_t sized  = 0;
    size_t sizes  = 0;
    size_t sizezc = 0;
    size_t sizeds = 0;

    /*
     * Get the number of element of each list.
     * Don't forget the -1 due to the last FUNCTION_NONE element
     */
#if defined(PRECISION_z)
    sizez  = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_z))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_c)
    sizec  = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_c))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_d)
    sized  = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_d))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_s)
    sizes  = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_s))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_zc)
    sizezc = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_zc)) / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_ds)
    sizeds = ( sizeof(PPTRACE_SYMBOL_LIST(coreblas_tcore_ds)) / sizeof(struct ezt_instrumented_function) ) - 1;
#endif

    /* Compute the total number of elements + the final FUNCTION_NONE element */
    size = sizez + sizec + sized + sizes + sizezc + sizeds + 1;

    PPTRACE_SYMBOL_LIST(coreblas_tcore) =
        malloc( sizeof( struct ezt_instrumented_function ) * size );

    /* Copy the lists with the global one */
    current_list_start = PPTRACE_SYMBOL_LIST(coreblas_tcore);

#if defined(PRECISION_z)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_z)),
            sizeof(struct ezt_instrumented_function) * sizez );
    current_list_start += sizez;
#endif

#if defined(PRECISION_c)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_c)),
            sizeof(struct ezt_instrumented_function) * sizec );
    current_list_start += sizec;
#endif

#if defined(PRECISION_d)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_d)),
            sizeof(struct ezt_instrumented_function) * sized );
    current_list_start += sized;
#endif

#if defined(PRECISION_s)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_s)),
            sizeof(struct ezt_instrumented_function) * sizes );
    current_list_start += sizes;
#endif

#if defined(PRECISION_zc)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_zc)),
            sizeof(struct ezt_instrumented_function) * sizezc );
    current_list_start += sizezc;
#endif

#if defined(PRECISION_ds)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(coreblas_tcore_ds)),
            sizeof(struct ezt_instrumented_function) * sizeds );
    current_list_start += sizeds;
#endif

    /* Add the final element */
    *current_list_start = (struct ezt_instrumented_function)FUNCTION_NONE_ELT;

    eztrace_log( dbg_lvl_debug, "eztrace_coreblas_tcore constructor starts\n" );
    EZT_REGISTER_MODULE( coreblas_tcore,
                         "Module for the coreblas_tcore library",
                         init_coreblas_tcore_,
                         finalize_coreblas_tcore );
    eztrace_log( dbg_lvl_debug, "eztrace_coreblas_tcore constructor ends\n" );
}

static void _coreblas_tcore_destroy( void ) __attribute__((destructor));
static void
_coreblas_tcore_destroy( void )
{
    free( PPTRACE_SYMBOL_LIST(coreblas_tcore) );
}
