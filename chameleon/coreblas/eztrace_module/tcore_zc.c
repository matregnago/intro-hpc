/**
 *
 * @file tcore_zc.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Eztrace module file generated from tcore_zc.h
 *
 * @version 1.4.0
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions mixed zc -> zc ds
 *
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>
#include "coreblas.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _coreblas_tcore_zc_initialized = 0;

#define CURRENT_MODULE coreblas_tcore_zc
DECLARE_CURRENT_MODULE;

/* pointers to actual coreblas_tcore_zc functions */

void (*libTCORE_clag2z)(int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) = NULL;
void (*libTCORE_zlag2c)(int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B, int *info) = NULL;



/* Wrapper functions */


void TCORE_clag2z(int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B) {
    FUNCTION_ENTRY;
    libTCORE_clag2z(M, N, A, B);
    FUNCTION_EXIT;

}

void TCORE_zlag2c(int M, int N, const CHAM_tile_t * A, CHAM_tile_t * B, int *info) {
    FUNCTION_ENTRY;
    libTCORE_zlag2c(M, N, A, B, info);
    FUNCTION_EXIT;

}


PPTRACE_START_INTERCEPT_FUNCTIONS(coreblas_tcore_zc)
INTERCEPT3("TCORE_clag2z", libTCORE_clag2z)
INTERCEPT3("TCORE_zlag2c", libTCORE_zlag2c)

PPTRACE_END_INTERCEPT_FUNCTIONS(coreblas_tcore_zc)
