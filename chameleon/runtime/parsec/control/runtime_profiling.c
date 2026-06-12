/**
 *
 * @file parsec/runtime_profiling.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC profiling routines
 *
 * @version 1.4.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @date 2024-02-18
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/timer.h"

double RUNTIME_get_time(){
    return CHAMELEON_timer();
}

void RUNTIME_start_profiling()
{
}

void RUNTIME_stop_profiling()
{
}

void RUNTIME_start_stats()
{
    chameleon_warning("RUNTIME_start_stats()", "pruning stats are not available with PaRSEC\n");
}

void RUNTIME_stop_stats()
{
    chameleon_warning("RUNTIME_stop_stats()", "pruning stats are not available with PaRSEC\n");
}

void RUNTIME_schedprofile_display(void)
{
    chameleon_warning("RUNTIME_schedprofile_display(parsec)", "Scheduler profiling is not available with PaRSEC\n");
}

void RUNTIME_kernelprofile_display(void)
{
    chameleon_warning("RUNTIME_kernelprofile_display(parsec)", "Kernel profiling is not available with PaRSEC\n");
}

/**
 *  Set iteration numbers for traces
 */
void RUNTIME_iteration_push( CHAM_context_t *chamctxt, unsigned long iteration )
{
    (void)chamctxt; (void)iteration;
    return;
}
void RUNTIME_iteration_pop( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

