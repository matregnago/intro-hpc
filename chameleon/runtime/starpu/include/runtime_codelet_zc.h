/**
 *
 * @file starpu/runtime_codelet_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU CHAMELEON_Complex64_t codelets header
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#ifndef _runtime_codelet_zc_h_
#define _runtime_codelet_zc_h_

#include <stdio.h>
#include "runtime_codelets.h"

#include "chameleon/tasks_zc.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_zc.h"
#include "coreblas/coreblas_zctile.h"
#if defined(CHAMELEON_USE_CUDA)
#include "gpucublas.h"
#endif
#endif

/*
 * MIXED PRECISION functions
 */
CODELETS_HEADER(zlag2c);
CODELETS_HEADER(clag2z);

#endif /* _runtime_codelet_zc_h_ */
