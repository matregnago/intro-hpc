/**
 *
 * @file global.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global coreblas variables and functions
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2025-12-19
 *
 */
#include "coreblas.h"
#include <stdarg.h>
#include <stdlib.h>
#include "chameleon/getenv.h"

int _coreblas_silent = 0;
static int coreblas_gemm3m_enabled = 0;

__attribute__((unused)) __attribute__((constructor)) static void
__coreblas_lib_init()
{
    _coreblas_silent = chameleon_getenv_get_value_int( "CHAMELEON_COREBLAS_SILENT", _coreblas_silent );
}

#if defined(CHAMELEON_KERNELS_TRACE)
void __coreblas_kernel_trace( const char *func, ... )
{
    char output[1024];
    int first = 1;
    int size = 0;
    int len = 1024;
    va_list va_list;
    const CHAM_tile_t *tile;

    if (_coreblas_silent) {
        return;
    }

    size += snprintf( output, len, "[coreblas] Execute %s(", func );

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

void
set_coreblas_gemm3m_enabled( int v ) {
    coreblas_gemm3m_enabled = v;
}

int
get_coreblas_gemm3m_enabled(void) {
    return coreblas_gemm3m_enabled;
}

/**
 *  LAPACK Constants
 */
char *chameleon_lapack_constants[] =
{
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 100
    "Row",                  // 101: ChamRowMajor
    "Column",               // 102: ChamColMajor
    "", "", "", "", "", "", "", "",
    "No transpose",         // 111: ChamNoTrans
    "Transpose",            // 112: ChamTrans
    "Conjugate transpose",  // 113: ChamConjTrans
    "", "", "", "", "", "", "",
    "Upper",                // 121: ChamUpper
    "Lower",                // 122: ChamLower
    "All",                  // 123: ChamUpperLower
    "", "", "", "", "", "", "",
    "Non-unit",             // 131: ChamNonUnit
    "Unit",                 // 132: ChamUnit
    "", "", "", "", "", "", "", "",
    "Left",                 // 141: ChamLeft
    "Right",                // 142: ChamRight
    "", "", "", "", "", "", "", "",
    "",                     // 151:
    "",                     // 152:
    "",                     // 153:
    "",                     // 154:
    "",                     // 155:
    "",                     // 156:
    "Epsilon",              // 157: ChamEps
    "",                     // 158:
    "",                     // 159:
    "",                     // 160:
    "", "", "", "", "", "", "", "", "", "",
    "One norm",             // 171: ChamOneNorm
    "",                     // 172: ChamRealOneNorm
    "",                     // 173: ChamTwoNorm
    "Frobenius norm",       // 174: ChamFrobeniusNorm
    "Infinity norm",        // 175: ChamInfNorm
    "",                     // 176: ChamRealInfNorm
    "Maximum norm",         // 177: ChamMaxNorm
    "",                     // 178: ChamRealMaxNorm
    "",                     // 179
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 200
    "Uniform",              // 201: ChamDistUniform
    "Symmetric",            // 202: ChamDistSymmetric
    "Normal",               // 203: ChamDistNormal
    "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 240
    "Hermitian",            // 241 ChamHermGeev
    "Positive ev Hermitian",// 242 ChamHermPoev
    "NonSymmetric pos sv",  // 243 ChamNonsymPosv
    "Symmetric pos sv",     // 244 ChamSymPosv
    "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 290
    "No Packing",           // 291 ChamNoPacking
    "U zero out subdiag",   // 292 ChamPackSubdiag
    "L zero out superdiag", // 293 ChamPackSupdiag
    "C",                    // 294 ChamPackColumn
    "R",                    // 295 ChamPackRow
    "B",                    // 296 ChamPackLowerBand
    "Q",                    // 297 ChamPackUpeprBand
    "Z",                    // 298 ChamPackAll
    "",                     // 299

    "",                     // 300
    "No vectors",           // 301 ChamNoVec
    "Vectors needed",       // 302 ChamVec, ChamSVDvrange
    "I",                    // 303 ChamIvec, ChamSVDirange
    "A",                    // 304 ChamAllVec, ChamSVDall
    "S",                    // 305 ChamSVec
    "O",                    // 306 ChamOVec
    "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 390
    "Forward",              // 391
    "Backward",             // 392
    "", "", "", "", "", "", "", "",
    "Columnwise",           // 401
    "Rowwise",              // 402
    "", "", "", "", "", "", "", ""  // Remember to add a coma!
};
