/**
 *
 * @file quark/chameleon_quark.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark runtime main header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_quark_h_
#define _chameleon_quark_h_

#include "control/common.h"

#include <quark.h>
#if defined(CHAMELEON_RUNTIME_SYNC)
#define QUARK_Insert_Task QUARK_Execute_Task
#endif
#include "coreblas.h"
#include "core_blas_dag.h"

typedef struct quark_option_s {
    Quark_Task_Flags flags;
    Quark *quark;
} quark_option_t;

static inline int cham_to_quark_access( cham_access_t accessA ) {
    if ( accessA == ChamR ) {
        return INPUT;
    }
    if ( accessA == ChamW ) {
        return OUTPUT;
    }
    return INOUT;
}

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( RUNTIME_data_getaddr( desc, m, n ) )

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _chameleon_quark_h_ */
