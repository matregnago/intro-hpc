/**
 *
 * @file coreblas.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU kernels main header
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#ifndef _coreblas_h_
#define _coreblas_h_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "coreblas/cblas_wrapper.h"

/**
 * CHAMELEON types and constants
 */
#include "chameleon/config.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/constants.h"

/**
 * CORE BLAS headers
 */
BEGIN_C_DECLS

#include "coreblas/coreblas_z.h"
#include "coreblas/coreblas_d.h"
#include "coreblas/coreblas_c.h"
#include "coreblas/coreblas_s.h"
#include "coreblas/coreblas_zc.h"
#include "coreblas/coreblas_ds.h"

#include "coreblas/coreblas_ztile.h"
#include "coreblas/coreblas_dtile.h"
#include "coreblas/coreblas_ctile.h"
#include "coreblas/coreblas_stile.h"
#include "coreblas/coreblas_zctile.h"
#include "coreblas/coreblas_dstile.h"

END_C_DECLS

/**
 * Coreblas Error
 */
#define coreblas_error(k, str) do {                                     \
        fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str) ;  \
        assert(0);                                                      \
    } while(0)

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

void set_coreblas_gemm3m_enabled( int v );
int  get_coreblas_gemm3m_enabled( void );


#if defined(CHAMELEON_KERNELS_TRACE)

void __coreblas_kernel_trace( const char *func, ... );
#define coreblas_kernel_trace( ... )                                    \
    do {                                                                \
        __coreblas_kernel_trace( __chameleon_func__, ##__VA_ARGS__, NULL );  \
    } while(0)

#else

#define coreblas_kernel_trace( ... ) do {} while(0)

#endif

void CORE_ipiv_to_perm( int m0, int m, int k, int K1, int K2, const int *ipiv, int *perm, int *invp );
void CORE_perm_to_idx( int m0, int m, int mt, const int *perm, const int *invp, int *perm_idx, int *invp_idx );

END_C_DECLS

#endif /* _coreblas_h_ */
