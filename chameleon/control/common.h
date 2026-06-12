/**
 *
 * @file common.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon common header file
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
/**
 *  CHAMELEON facilities of interest to both CHAMELEON core developer
 *  and also of interest to CHAMELEON community contributor.
 */
#ifndef _chameleon_common_h_
#define _chameleon_common_h_

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE 1
#endif
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/**
 *  Chameleon header files
 */
#include "chameleon.h"

#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

/**
 * Implementation headers
 */
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#if defined(CHAMELEON_USE_OPENCL) && !defined(CHAMELEON_SIMULATION)
#include <OpenCL/cl.h>
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

#include "control/global.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "control/descriptor.h"
#include "control/async.h"
#include "chameleon/getenv.h"

/**
 *  Global shortcuts
 */
#define CHAMELEON_RANK        chameleon_rank(chamctxt)
#define CHAMELEON_NB          chamctxt->nb
#define CHAMELEON_IB          chamctxt->ib
#define CHAMELEON_RHBLK       chamctxt->rhblock
#define CHAMELEON_TRANSLATION chamctxt->translation
#define CHAMELEON_PARALLEL    chamctxt->parallel_enabled
#define CHAMELEON_STATISTICS  chamctxt->statistics_enabled

/**
 *  IPT internal define
 */
#define ChamIPT_NoDep   0
#define ChamIPT_Panel   1
#define ChamIPT_All     2

/**
 * Structure for reduction operations
 */
struct chameleon_reduce_s {
    cham_getrf_allreduce_t  alg_allreduce; /**< Specifies the algorithm used for the allreduce                          */
    int                    *proc_involved; /**< Specifies the processes involved in the reduction operation             */
    unsigned int            involved;      /**< Specifies if the current process is involved in the reduction operation */
    int                     np_involved;   /**< Specifies the number of involved processes in the reduction operation   */
    int                     arity;         /**< Specifies the arity of the reduction tree                               */
};
typedef struct chameleon_reduce_s CHAM_reduce_t;

/**
 *  Global array of LAPACK constants
 */
extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GNUC__)
static inline int chameleon_asprintf( char **restrict strp, const char *fmt, ... ) __attribute__((format(printf,2,3)));
#endif
static inline int chameleon_asprintf( char **restrict strp, const char *fmt, ... )
{
    va_list ap;
    int rc;

    va_start( ap, fmt );
    rc = vasprintf( strp, fmt, ap );
    va_end( ap );

    assert( rc != -1 );
    return rc;
}

#if defined(CHAMELEON_KERNELS_TRACE)
static inline const char *chameleon_codelet_name( const char *kname, int nbtiles, ... )
{
    int  nbchar = 256;
    char output[256];
    char *strptr = output;
    va_list ap;
    int rc, i;

    rc = snprintf( strptr, nbchar, "%s(", kname );
    if ( rc < 0 ) {
        return NULL;
    }
    strptr += rc;
    nbchar -= rc;

    va_start( ap, nbtiles );
    for( i=0; i<nbtiles; i++ )
    {
        CHAM_tile_t *tile = va_arg( ap, CHAM_tile_t* );
        if ( i > 0 ) {
            strptr[0] = ',';
            strptr++;
            nbchar--;
        }
        rc = snprintf( strptr, nbchar, " %s", tile->name );
        if ( rc < 0 ) {
            return NULL;
        }
        strptr += rc;
        nbchar -= rc;
    }
    va_end( ap );

    if ( nbchar < 3 ) {
        return NULL;
    }
    strptr[0] = ' ';
    strptr[1] = ')';
    strptr[2] = '\0';

    return strndup( output, 255 );
}
#else
static inline const char *chameleon_codelet_name( const char *kname, int nbtiles, ... )
{
    (void)nbtiles;
    return kname;
}
#endif

void chameleon_pmap( cham_uplo_t uplo, int ndata, cham_map_data_t *data,
                     cham_map_operator_t *op_fct, void *op_args,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

#include "control/compute_s.h"
#include "control/compute_d.h"
#include "control/compute_c.h"
#include "control/compute_z.h"

/*
void chameleon_pdlag2s(CHAM_context_t *chamctxt);
void chameleon_pzlag2c(CHAM_context_t *chamctxt);
void chameleon_pslag2d(CHAM_context_t *chamctxt);
void chameleon_pclag2z(CHAM_context_t *chamctxt);
*/

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_common_h_ */
