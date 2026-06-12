/**
 *
 * @file cudaglobal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global gpucublas variables and functions
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#include "gpucublas.h"
#include <stdarg.h>
#include <stdlib.h>
#include "chameleon/getenv.h"

int _gpucublas_silent = 0;

__attribute__((unused)) __attribute__((constructor)) static void
__gpucublas_lib_init()
{
    _gpucublas_silent = chameleon_getenv_get_value_int( "CHAMELEON_GPUCUBLAS_SILENT", _gpucublas_silent );
}

#if defined(CHAMELEON_KERNELS_TRACE)
void __gpucublas_kernel_trace( const char *func, ... )
{
    char output[1024];
    int first = 1;
    int size = 0;
    int len = 1024;
    va_list va_list;
    const CHAM_tile_t *tile;

    if (_gpucublas_silent) {
        return;
    }

    size += snprintf( output, len, "[gpucublas] Execute %s(", func );

    va_start( va_list, func );
    while((tile = va_arg(va_list, const CHAM_tile_t*)) != 0) {
        size += snprintf( output+size, len-size, "%s%s",
                          first ? "" : ", ",
                          tile->name );
        size += snprintf( output+size, len-size, " / %p",
                          CHAM_tile_get_ptr( tile ) );
        first = 0;
    }
    va_end( va_list );

    fprintf( stderr, "%s)\n", output );
    fflush(stderr);
}
#endif

/**
 *  LAPACK Constants
 */
int chameleon_cublas_constants[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 100
    0,                      // 101: ChamRowMajor
    0,                      // 102: ChamColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: ChamNoTrans
    CUBLAS_OP_T,            // 112: ChamTrans
    CUBLAS_OP_C,            // 113: ChamConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: ChamUpper
    CUBLAS_FILL_MODE_LOWER, // 122: ChamLower
    CUBLAS_FILL_MODE_FULL,  // 123: ChamUpperLower
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: ChamNonUnit
    CUBLAS_DIAG_UNIT,       // 132: ChamUnit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: ChamLeft
    CUBLAS_SIDE_RIGHT,      // 142: ChamRight
    0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 151:
    0,                      // 152:
    0,                      // 153:
    0,                      // 154:
    0,                      // 155:
    0,                      // 156:
    0,                      // 157: ChamEps
    0,                      // 158:
    0,                      // 159:
    0,                      // 160:
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 171: ChamOneNorm
    0,                      // 172: ChamRealOneNorm
    0,                      // 173: ChamTwoNorm
    0,                      // 174: ChamFrobeniusNorm
    0,                      // 175: ChamInfNorm
    0,                      // 176: ChamRealInfNorm
    0,                      // 177: ChamMaxNorm
    0,                      // 178: ChamRealMaxNorm
    0,                      // 179
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 200
    0,                      // 201: ChamDistUniform
    0,                      // 202: ChamDistSymmetric
    0,                      // 203: ChamDistNormal
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 240
    0,                      // 241 ChamHermGeev
    0,                      // 242 ChamHermPoev
    0,                      // 243 ChamNonsymPosv
    0,                      // 244 ChamSymPosv
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 290
    0,                      // 291 ChamNoPacking
    0,                      // 292 ChamPackSubdiag
    0,                      // 293 ChamPackSupdiag
    0,                      // 294 ChamPackColumn
    0,                      // 295 ChamPackRow
    0,                      // 296 ChamPackLowerBand
    0,                      // 297 ChamPackUpeprBand
    0,                      // 298 ChamPackAll
    0,                      // 299
    0,                      // 300
    0,                      // 301 ChamNoVec
    0,                      // 302 ChamVec, ChamSVDvrange
    0,                      // 303 ChamIvec, ChamSVDirange
    0,                      // 304 ChamAllVec, ChamSVDall
    0,                      // 305 ChamSVec
    0,                      // 306 ChamOVec
    0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 390
    0,                      // 391 Forward
    0,                      // 392 Backward
    0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 401 Columnwise
    0,                      // 402 Rowwise
    0, 0, 0, 0, 0, 0, 0, 0  // Remember to add a coma!
};
