/**
 *
 * @file starpu/codelet_convert.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-09-17
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelets.h"
#define PRECISION_z
#include "runtime_codelet_z.h"
#undef PRECISION_z
#define PRECISION_d
#include "runtime_codelet_d.h"
#undef PRECISION_d
#define PRECISION_c
#include "runtime_codelet_c.h"
#undef PRECISION_c
#define PRECISION_s
#include "runtime_codelet_s.h"
#undef PRECISION_s
#define PRECISION_zc
#include "runtime_codelet_zc.h"
#undef PRECISION_zc
#define PRECISION_ds
#include "runtime_codelet_ds.h"
#undef PRECISION_ds

void
insert_task_convert( const RUNTIME_option_t *options,
                     int m, int n,
                     cham_flttype_t       fromtype,
                     starpu_data_handle_t fromtile,
                     cham_flttype_t       totype,
                     starpu_data_handle_t totile )
{
    struct starpu_codelet *codelet = NULL;
    callback_fct_t cl_convert_callback = NULL;

    int conversion = ChamConvert( fromtype, totype );

    switch( conversion ) {
#if defined(CHAMELEON_PREC_ZC)
    case ChamConvertComplexDoubleToSingle:
        codelet = &cl_zlag2c;
        cl_convert_callback = cl_zlag2c_callback;
        break;

    case ChamConvertComplexSingleToDouble:
        codelet = &cl_clag2z;
        cl_convert_callback = cl_clag2z_callback;
        break;
#endif

#if defined(CHAMELEON_PREC_DS)
    case ChamConvertRealDoubleToSingle:
        codelet = &cl_dlag2s;
        cl_convert_callback = cl_dlag2s_callback;
        break;

    case ChamConvertRealSingleToDouble:
        codelet = &cl_slag2d;
        cl_convert_callback = cl_slag2d_callback;
        break;
#endif

#if defined(CHAMELEON_PREC_D) && defined(GPUCUBLAS_HAVE_CUDA_HALF)
    case ChamConvertRealDoubleToHalf:
        codelet = &cl_dlag2h;
        cl_convert_callback = cl_dlag2h_callback;
        break;

    case ChamConvertRealHalfToDouble:
        codelet = &cl_hlag2d;
        cl_convert_callback = cl_hlag2d_callback;
        break;
#endif

#if defined(CHAMELEON_PREC_S) && defined(GPUCUBLAS_HAVE_CUDA_HALF)
    case ChamConvertRealSingleToHalf:
        codelet = &cl_slag2h;
        cl_convert_callback = cl_slag2h_callback;
        break;

    case ChamConvertRealHalfToSingle:
        codelet = &cl_hlag2s;
        cl_convert_callback = cl_hlag2s_callback;
        break;
#endif

    case ChamConvertComplexDoubleToDouble:
        return;
    case ChamConvertComplexSingleToSingle:
        return;
    case ChamConvertComplexHalfToHalf:
        return;
    case ChamConvertRealDoubleToDouble:
        return;
    case ChamConvertRealSingleToSingle:
        return;
    case ChamConvertRealHalfToHalf:
        return;

    default:
        assert(0);
        fprintf( stderr, "INSERT_TASK_convert: Unknown conversion type\n" );
        return;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        codelet,
        /* Task codelet arguments */
        STARPU_VALUE, &m, sizeof(int),
        STARPU_VALUE, &n, sizeof(int),

        /* Task handles */
        STARPU_R, fromtile,
        STARPU_W, totile,

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( convert ),
        0);

    return;
}
