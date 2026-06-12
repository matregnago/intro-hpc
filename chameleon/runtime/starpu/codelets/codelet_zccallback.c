/**
 *
 * @file starpu/codelet_zccallback.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zc callback StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Florent Pruvost
 * @date 2024-09-17
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_zc.h"

CHAMELEON_CL_CB(zlag2c,        cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), 0,                                      M*N)
CHAMELEON_CL_CB(clag2z,        cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), 0,                                      M*N)
