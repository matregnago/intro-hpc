/**
 *
 * @file chameleon_simulate.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon SimGrid simulation header
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_simulate_h_
#define _chameleon_simulate_h_

#include "chameleon/config.h"

/* we need this when starpu is compiled with simgrid enabled */
#if defined(CHAMELEON_SIMULATION)

#if !defined(CHAMELEON_SCHED_STARPU)
#error "CHAMELEON_SIMULATION can not be enabled without StarPU support"
#endif

#ifndef main
#define main starpu_main
#endif

#if defined(CHAMELEON_USE_MPI)
#define starpu_main smpi_simulated_main_
#endif

#endif

#endif /* _chameleon_simulate_h_ */
