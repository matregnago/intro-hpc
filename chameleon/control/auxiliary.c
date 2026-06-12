/**
 *
 * @file auxiliary.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 ***
 *
 * @defgroup Auxiliary
 * @brief Group auxiliary routines exposed to users
 *
 */

#include "control/common.h"
#include "control/auxiliary.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int _chameleon_silent = 0;

__attribute__((unused)) __attribute__((constructor)) static void
__chameleon_lib_init()
{
    _chameleon_silent = chameleon_getenv_get_value_int( "CHAMELEON_SILENT", _chameleon_silent );
}

void
__chameleon_log_debug( const char *fmt, ... )
{
    va_list ap;

    if ( _chameleon_silent ) {
        return;
    }

    va_start( ap, fmt );
    vfprintf( stderr, fmt, ap );
    va_end(ap);
}

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while CHAMELEON is being used correctly.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void chameleon_warning(const char *func_name, const char *msg_text)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("chameleon_warning", "CHAMELEON not initialized");
    }
    if (chamctxt->warnings_enabled) {
        fprintf(stderr, "CHAMELEON WARNING: %s(): %s\n", func_name, msg_text);
    }
}

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of CHAMELEON.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void chameleon_error(const char *func_name, const char *msg_text)
{
    fprintf( stderr, "CHAMELEON ERROR: %s(): %s\n", func_name, msg_text );
}

/**
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *  Context oblivious.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void chameleon_fatal_error(const char *func_name, const char *msg_text)
{
    fprintf(stderr, "CHAMELEON FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/**
 *  Returns core id
 */
int chameleon_rank(CHAM_context_t *chamctxt)
{
    return RUNTIME_thread_rank( chamctxt );
}

/**
 *  Tune block size nb and internal block size ib
 */
int chameleon_tune(cham_tasktype_t func, int M, int N, int NRHS)
{
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if ( chamctxt && chamctxt->autotuning_enabled == CHAMELEON_TRUE ) {
        chameleon_warning( "chameleon_tune", "Autotunning not available for now" );
    }
    (void)func;
    (void)M;
    (void)N;
    (void)NRHS;
    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  CHAMELEON_Version - Reports CHAMELEON version number.
 *
 ******************************************************************************
 *
 * @param[out] ver_major
 *          CHAMELEON major version number.
 *
 * @param[out] ver_minor
 *          CHAMELEON minor version number.
 *
 * @param[out] ver_micro
 *          CHAMELEON micro version number.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Version(int *ver_major, int *ver_minor, int *ver_micro)
{
    if (! ver_major && ! ver_minor && ! ver_micro)
        return  CHAMELEON_ERR_ILLEGAL_VALUE;

    if (ver_major)
        *ver_major = CHAMELEON_VERSION_MAJOR;

    if (ver_minor)
        *ver_minor = CHAMELEON_VERSION_MINOR;

    if (ver_micro)
        *ver_micro = CHAMELEON_VERSION_MICRO;

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  CHAMELEON_Element_Size - Reports the size in bytes of a CHAMELEON precision type
 *  (e.g. ChamInteger, ChamRealFloat, etc).
 *
 ******************************************************************************
 *
 * @param[in] type
 *          CHAMELEON element type, can be one of the following:
 *          - ChamByte
 *          - ChamInteger
 *          - ChamRealFloat
 *          - ChamRealDouble
 *          - ChamComplexFloat
 *          - ChamComplexDouble
 *
 ******************************************************************************
 *
 * @retval Element size in bytes
 *
 */
ssize_t CHAMELEON_Element_Size( cham_flttype_t type )
{
    switch( cham_get_flttype(type) ) {
        case ChamByte:          return (ssize_t)1;
        case ChamInteger16:     return sizeof(int16_t);
        case ChamInteger32:     return sizeof(int32_t);
        case ChamInteger64:     return sizeof(int64_t);
        case ChamRealHalf:      return (ssize_t)2;
        case ChamRealFloat:     return sizeof(float);
        case ChamRealDouble:    return sizeof(double);
        case ChamComplexHalf:   return (ssize_t)4;
        case ChamComplexFloat:  return sizeof(float)  * 2;
        case ChamComplexDouble: return sizeof(double) * 2;
        default: chameleon_fatal_error("CHAMELEON_Element_Size", "undefined type");
                 return CHAMELEON_ERR_ILLEGAL_VALUE;

    }
}

/**
 *  Display a progress percentage in stderr
 */
void update_progress(int currentValue, int maximumValue) {
    div_t res ;
    static int progress = -1; /* varie de 0 a 100 au cours du calcul concerne */

    if (maximumValue == 0) {
        res.quot = 100;
    }
    else {
        if (currentValue < (INT_MAX / 100) ) {
            res = div(currentValue*100, maximumValue);
        }
        else {
            /* Compute the quotient */
            res.quot = (int)( ((double)currentValue * 100.) / (double)maximumValue );
        }
    }

    // Print the percentage
    if (res.quot > progress) {
        fprintf(stderr, "%3d%%\b\b\b\b", res.quot);
    }
    progress = res.quot;

    if (currentValue >= maximumValue) {
        progress = -1;
    }
}

// A function to display the progress indicator.
// By default it is update_progress()
// The user can change it with CHAMELEON_Set_update_progress_callback()
void (*update_progress_callback)(int, int) = update_progress;

int CHAMELEON_Set_update_progress_callback(void (*p)(int, int)) {
  update_progress_callback = p;
  return CHAMELEON_SUCCESS;
}

