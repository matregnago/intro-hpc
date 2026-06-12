/**
 *
 * @file context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context management routines
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Alycia Lisito
 * @author Matthieu Kuhn
 * @author Loris Lucido
 * @author Terry Cojean
 * @author Xavier Lacoste
 * @date 2025-12-19
 *
 ***
 *
 * @defgroup Options
 * @brief Group routines exposed to users to handle options
 *
 */
#include "control/common.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "chameleon/runtime.h"
#include <stdlib.h>
#include <limits.h>

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas.h"
#endif

/**
 *  Global data
 */
/* master threads context lookup table */
static CHAM_context_t *chameleon_ctxt = NULL;

static inline cham_householder_t
chameleon_getenv_householder(char * string, cham_householder_t default_value) {
    long int ret;
    char *str = chameleon_getenv(string);
    if (str == NULL) return default_value;

    if ( sscanf( str, "%ld", &ret ) == 1 ) {
        switch (ret) {
        case ChamFlatHouseholder:
            return ChamFlatHouseholder;
        case ChamTreeHouseholder:
            return ChamTreeHouseholder;
        default:
            chameleon_error( "chameleon_getenv_householder", "Incorrect householder value" );
            return default_value;
        }
    }

    if(0 == strcasecmp("chamflathouseholder", str)) { return ChamFlatHouseholder; }
    if(0 == strcasecmp("flat", str)) { return ChamFlatHouseholder; }
    if(0 == strcasecmp("chamtreehouseholder", str)) { return ChamTreeHouseholder; }
    if(0 == strcasecmp("tree", str)) { return ChamTreeHouseholder; }

    chameleon_error( "chameleon_getenv_householder", "Incorrect householder type" );

    return default_value;
}

static inline cham_translation_t
chameleon_getenv_translation(char * string, cham_translation_t default_value) {
    long int ret;
    char *str = chameleon_getenv(string);
    if (str == NULL) return default_value;

    if ( sscanf( str, "%ld", &ret ) == 1 ) {
        switch (ret) {
        case ChamInPlace:
            return ChamInPlace;
        case ChamOutOfPlace:
            return ChamOutOfPlace;
        default:
            chameleon_error( "chameleon_getenv_translation", "Incorrect translation value" );
            return default_value;
        }
    }

    if(0 == strcasecmp("chaminplace", str)) { return ChamInPlace; }
    if(0 == strcasecmp("inplace", str)) { return ChamInPlace; }
    if(0 == strcasecmp("in", str)) { return ChamInPlace; }

    if(0 == strcasecmp("chamoutofplace", str)) { return ChamOutOfPlace; }
    if(0 == strcasecmp("outofplace", str)) { return ChamOutOfPlace; }
    if(0 == strcasecmp("out", str)) { return ChamOutOfPlace; }

    chameleon_error( "chameleon_getenv_translation", "Incorrect translation type" );

    return default_value;
}

/**
 *  Create new context
 */
CHAM_context_t *chameleon_context_create()
{
    CHAM_context_t *chamctxt;

    if ( chameleon_ctxt != NULL ) {
        chameleon_error("chameleon_context_create", "a context is already existing\n");
        return NULL;
    }

    chamctxt = (CHAM_context_t*)malloc(sizeof(CHAM_context_t));
    if (chamctxt == NULL) {
        chameleon_error("chameleon_context_create", "malloc() failed");
        return NULL;
    }

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    chamctxt->nb                  = chameleon_getenv_get_value_int( "CHAMELEON_TILE_SIZE",        384 );
    chamctxt->ib                  = chameleon_getenv_get_value_int( "CHAMELEON_INNER_BLOCK_SIZE",  48 );
    chamctxt->rhblock             = chameleon_getenv_get_value_int( "CHAMELEON_HOUSEHOLDER_SIZE",   4 );
    chamctxt->lookahead           = chameleon_getenv_get_value_int( "CHAMELEON_LOOKAHEAD",          1 );

    chamctxt->nworkers            = 1;
    chamctxt->ncudas              = 0;
    chamctxt->nthreads_per_worker = 1;

    chamctxt->warnings_enabled    = chameleon_env_on_off( "CHAMELEON_WARNINGS",        CHAMELEON_TRUE  );
    chamctxt->autotuning_enabled  = chameleon_env_on_off( "CHAMELEON_AUTOTUNING",      CHAMELEON_FALSE );
    chamctxt->parallel_enabled    = chameleon_env_on_off( "CHAMELEON_PARALLEL_KERNEL", CHAMELEON_FALSE );
    chamctxt->statistics_enabled  = chameleon_env_on_off( "CHAMELEON_GENERATE_STATS",  CHAMELEON_FALSE );
    chamctxt->progress_enabled    = chameleon_env_on_off( "CHAMELEON_PROGRESS",        CHAMELEON_FALSE );
    chamctxt->generic_enabled     = chameleon_env_on_off( "CHAMELEON_GENERIC",         CHAMELEON_FALSE );
    chamctxt->autominmax_enabled  = chameleon_env_on_off( "CHAMELEON_AUTOMINMAX",      CHAMELEON_TRUE  );
    chamctxt->optlacpy_enabled    = chameleon_env_on_off( "CHAMELEON_OPTIMIZED_LACPY", CHAMELEON_TRUE  );

    chamctxt->runtime_paused      = CHAMELEON_FALSE;

    chamctxt->householder = chameleon_getenv_householder( "CHAMELEON_HOUSEHOLDER_MODE", ChamFlatHouseholder );
    chamctxt->translation = chameleon_getenv_translation( "CHAMELEON_TRANSLATION_MODE", ChamInPlace );

    /* First and last step for main algorithms to get partial traces for example */
    chamctxt->first_step = chameleon_getenv_get_value_int( "CHAMELEON_FIRST_STEP", 0       );
    chamctxt->first_step = chameleon_max( chamctxt->first_step, 0 );
    chamctxt->last_step  = chameleon_getenv_get_value_int( "CHAMELEON_LAST_STEP",  INT_MAX );
    chamctxt->last_step  = ( chamctxt->last_step >= 0 ) ? chamctxt->last_step : 0;

    /* Initialize scheduler */
    RUNTIME_context_create(chamctxt);

    chameleon_ctxt = chamctxt;
    return chamctxt;
}


/**
 *  Return context for a thread
 */
CHAM_context_t *chameleon_context_self()
{
    return chameleon_ctxt;
}

/**
 *  Clean the context
 */
int chameleon_context_destroy(){

    RUNTIME_context_destroy(chameleon_ctxt);
    free(chameleon_ctxt);
    chameleon_ctxt = NULL;

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Enable - Enable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be enabled:
 *          @arg CHAMELEON_WARNINGS       Printing of warning messages.
 *          @arg CHAMELEON_AUTOTUNING     Autotuning for tile size and inner block size.
 *          @arg CHAMELEON_GENERATE_TRACE Enable/start the trace generation.
 *          @arg CHAMELEON_GENERATE_STATS Enable/start the kernel statistics.
 *          @arg CHAMELEON_PROGRESS       Enable the progress indicator.
 *          @arg CHAMELEON_GEMM3M         Use z/cgemm3m for complexe matrix-matrix products.
 *          @arg CHAMELEON_GENERIC        Enforce generic algorithms instead of 2DBC specific algotihms.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Enable(int option)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Enable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case CHAMELEON_WARNINGS:
            chamctxt->warnings_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_AUTOTUNING:
            chamctxt->autotuning_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_GENERATE_TRACE:
            RUNTIME_start_profiling();
            break;
        case CHAMELEON_GENERATE_STATS:
            chamctxt->statistics_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_PROGRESS:
            chamctxt->progress_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(1);
#else
            chameleon_error("CHAMELEON_Enable", "cannot enable GEMM3M (not available in cblas)");
#endif
            break;
        case CHAMELEON_PARALLEL_KERNEL:
            chamctxt->parallel_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_GENERIC:
            chamctxt->generic_enabled = CHAMELEON_TRUE;
            break;
        default:
            chameleon_error("CHAMELEON_Enable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
        case CHAMELEON_BOUND:
            break;
    }

    /* Enable at the lower level if required */
    RUNTIME_enable( chamctxt->schedopt, option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Disable - Disable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be disabled:
 *          @arg CHAMELEON_WARNINGS       Printing of warning messages.
 *          @arg CHAMELEON_AUTOTUNING     Autotuning for tile size and inner block size.
 *          @arg CHAMELEON_GENERATE_TRACE Disable/pause the trace generation.
 *          @arg CHAMELEON_GENERATE_STATS Disable/pause the kernel statistics.
 *          @arg CHAMELEON_PROGRESS       Disable the progress indicator.
 *          @arg CHAMELEON_GEMM3M         Disabel the use of z/cgemm3m for complexe matrix-matrix products.
 *          @arg CHAMELEON_GENERIC        Do not enforce generic algorithms instead of 2DBC specific algotihms.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Disable(int option)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Disable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case CHAMELEON_WARNINGS:
            chamctxt->warnings_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_AUTOTUNING:
            chamctxt->autotuning_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GENERATE_TRACE:
            RUNTIME_stop_profiling();
            break;
        case CHAMELEON_GENERATE_STATS:
            chamctxt->statistics_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_PROGRESS:
            chamctxt->progress_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(0);
#endif
            break;
        case CHAMELEON_PARALLEL_KERNEL:
            chamctxt->parallel_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GENERIC:
            chamctxt->generic_enabled = CHAMELEON_FALSE;
            break;
        default:
            chameleon_error("CHAMELEON_Disable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Disable at the lower level if required */
    RUNTIME_disable( chamctxt->schedopt, option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Set - Set CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg CHAMELEON_LOOKAHEAD:        depth of the look ahead in algorithms
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Set( int param, int value )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Set", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative tile size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->nb = value;
            if ( chamctxt->autotuning_enabled ) {
                chamctxt->autotuning_enabled = CHAMELEON_FALSE;
                chameleon_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            chamctxt->ib = chameleon_min( chamctxt->nb, chamctxt->ib );
            break;
        case CHAMELEON_INNER_BLOCK_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative inner block size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            if (value > chamctxt->nb) {
                chameleon_warning("CHAMELEON_Set", "inner block larger than tile");
            }
            chamctxt->ib = value;
            /* Limit ib to nb */
            chamctxt->ib = chameleon_min( chamctxt->nb, chamctxt->ib );

            if ( chamctxt->autotuning_enabled ) {
                chamctxt->autotuning_enabled = CHAMELEON_FALSE;
                chameleon_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            break;
        case CHAMELEON_HOUSEHOLDER_MODE:
            if (value != ChamFlatHouseholder && value != ChamTreeHouseholder) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_HOUSEHOLDER_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->householder = value;
            break;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative householder size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->rhblock = value;
            break;
        case CHAMELEON_TRANSLATION_MODE:
            if (value != ChamInPlace && value != ChamOutOfPlace) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_TRANSLATION_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->translation = value;
            break;
        case CHAMELEON_LOOKAHEAD:
            if (value < 1) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_LOOKAHEAD");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->lookahead = value;
            break;
        default:
            chameleon_error("CHAMELEON_Set", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Get - Get value of CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg CHAMELEON_LOOKAHEAD:        depth of the look ahead in algorithms,
 *          @arg CHAMELEON_HOUSEHOLDER_MODE: ChamFlatHouseholder or ChamTreeHouseholder,
 *          @arg CHAMELEON_HOUSEHOLDER_SIZE: Size of the Householder tree,
 *          @arg CHAMELEON_TRANSLATION_MODE: ChamInPlace or ChamOutOfPlace,
 *          @arg CHAMELEON_LOOKAHEAD:        Size of the look-ahead,
 *          @arg CHAMELEON_RUNTIME:          Runtime used by the library (See RUNTIME_id_t)
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Get( int param, int *value )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Get", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            *value = chamctxt->nb;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_INNER_BLOCK_SIZE:
            *value = chamctxt->ib;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_MODE:
            *value = chamctxt->householder;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            *value = chamctxt->rhblock;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_TRANSLATION_MODE:
            *value = chamctxt->translation;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_LOOKAHEAD:
            *value = chamctxt->lookahead;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_RUNTIME:
            *value = chamctxt->scheduler;
            return CHAMELEON_SUCCESS;
        default:
            chameleon_error("CHAMELEON_Get", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_GetContext - Get the Chameleon context pointer
 *
 *******************************************************************************
 *
 * @retval NULL if not initialized
 * @retval The pointer to the initialized Chameleon context
 *
 */
CHAM_context_t *CHAMELEON_GetContext()
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Get", "CHAMELEON not initialized");
        return NULL;
    }
    return chamctxt;
}
