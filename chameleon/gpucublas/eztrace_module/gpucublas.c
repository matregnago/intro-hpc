/**
 *
 * @file gpucublas.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon gpucublas EZtrace 2.0 module
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
#define FUNCTION_NONE_ELT  {			\
        .function_name = "",                    \
        .callback      = NULL,                  \
        .event_id      = -1,                    \
    }

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _gpucublas_initialized = 0;

struct ezt_instrumented_function *PPTRACE_SYMBOL_LIST(gpucublas);

#if defined(PRECISION_z)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_z);
#endif
#if defined(PRECISION_c)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_c);
#endif
#if defined(PRECISION_d)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_d);
#endif
#if defined(PRECISION_s)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_s);
#endif

#if defined(PRECISION_zc)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_zc);
#endif
#if defined(PRECISION_ds)
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(gpucublas_ds);
#endif

static void init_gpucublas_()
{
    INSTRUMENT_FUNCTIONS(gpucublas);

    if (eztrace_autostart_enabled())
        eztrace_start();

    _gpucublas_initialized = 1;
}

static void finalize_gpucublas()
{
    _gpucublas_initialized = 0;

    eztrace_stop();
}

static void _gpucublas_init(void) __attribute__((constructor));
static void _gpucublas_init(void)
{
    struct ezt_instrumented_function *current_list_start;
    size_t size   = 0;
    size_t sizez  = 0;
    size_t sizec  = 0;
    size_t sized  = 0;
    size_t sizes  = 0;
    size_t sizezc = 0;
    size_t sizeds = 0;
    size_t sizeh  = 0;

    /*
     * Get the number of element of each list.
     * Don't forget the -1 due to the last FUNCTION_NONE element
     */
#if defined(PRECISION_z)
    sizez  = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_z))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_c)
    sizec  = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_c))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_d)
    sized  = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_d))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_s)
    sizes  = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_s))  / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
#if defined(PRECISION_zc)
    sizezc = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_zc)) / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#if defined(PRECISION_ds)
    sizeds = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_ds)) / sizeof(struct ezt_instrumented_function) ) - 1;
#endif
#endif

#if defined(GPUCUBLAS_HAVE_CUBLASHGEMM) || defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
    sizeh = ( sizeof(PPTRACE_SYMBOL_LIST(gpucublas_h))   / sizeof(struct ezt_instrumented_function) ) - 1;
#endif

    /* Compute the total number of elements + the final FUNCTION_NONE element */
    size = sizez + sizec + sized + sizes + sizezc + sizeds + sizeh + 1;

    PPTRACE_SYMBOL_LIST(gpucublas) = malloc( sizeof(struct ezt_instrumented_function) * size );

    /* Copy the lists with the global one */
    current_list_start = PPTRACE_SYMBOL_LIST(gpucublas);

#if defined(PRECISION_z)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_z)),
            sizeof(struct ezt_instrumented_function) * sizez );
    current_list_start += sizez;
#endif

#if defined(PRECISION_c)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_c)),
            sizeof(struct ezt_instrumented_function) * sizec );
    current_list_start += sizec;
#endif

#if defined(PRECISION_d)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_d)),
            sizeof(struct ezt_instrumented_function) * sized );
    current_list_start += sized;
#endif

#if defined(PRECISION_s)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_s)),
            sizeof(struct ezt_instrumented_function) * sizes );
    current_list_start += sizes;
#endif

#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
#if defined(PRECISION_zc)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_zc)),
            sizeof(struct ezt_instrumented_function) * sizezc );
    current_list_start += sizezc;
#endif

#if defined(PRECISION_ds)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_ds)),
            sizeof(struct ezt_instrumented_function) * sizeds );
    current_list_start += sizeds;
#endif
#endif

#if defined(GPUCUBLAS_HAVE_CUBLASHGEMM) || defined(GPUCUBLAS_HAVE_CUBLASGEMMEX)
    memcpy( (void *)current_list_start,
            (void *)&(PPTRACE_SYMBOL_LIST(gpucublas_h)),
            sizeof(struct ezt_instrumented_function) * sizeh );
    current_list_start += sizeh;
#endif

    /* Add the final element */
    *current_list_start = (struct ezt_instrumented_function)FUNCTION_NONE_ELT;

    eztrace_log(dbg_lvl_debug, "eztrace_gpucublas constructor starts\n");
    EZT_REGISTER_MODULE(gpucublas, "Module for the gpucublas library",
                        init_gpucublas_, finalize_gpucublas);
    eztrace_log(dbg_lvl_debug, "eztrace_gpucublas constructor ends\n");
}

static void _gpucublas_destroy(void) __attribute__((destructor));
static void _gpucublas_destroy(void)
{
    free(PPTRACE_SYMBOL_LIST(gpucublas));
}
