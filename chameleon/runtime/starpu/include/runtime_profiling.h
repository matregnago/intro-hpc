/**
 *
 * @file starpu/runtime_profiling.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU profiling and kernel locality header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-12-19
 *
 */
#ifndef _runtime_profiling_h_
#define _runtime_profiling_h_

#ifdef CHAMELEON_ENABLE_PRUNING_STATS
extern unsigned long RUNTIME_total_tasks;
extern unsigned long RUNTIME_exec_tasks;
extern unsigned long RUNTIME_comm_tasks;
extern unsigned long RUNTIME_changed_tasks;
#endif

typedef struct measure_s {
    double sum;
    double sum2;
    long   n;
} measure_t;

void RUNTIME_profiling_display_info(const char *kernel_name, measure_t perf[STARPU_NMAXWORKERS]);
void RUNTIME_profiling_display_efficiency(void);

void RUNTIME_profiling_zdisplay_all(void);
void RUNTIME_profiling_cdisplay_all(void);
void RUNTIME_profiling_ddisplay_all(void);
void RUNTIME_profiling_sdisplay_all(void);

#endif /* _runtime_profiling_h_ */
