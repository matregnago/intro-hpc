/**
 *
 * @file chameleon_types.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon basic datatypes header
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Lucas Barros de Assis
 * @author Thomas Mijieux
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_types_h_
#define _chameleon_types_h_

#include "chameleon/config.h"
#include <math.h>

/**
 *  System requirements
 */
#include <stddef.h>
#if defined( _WIN32 )
  /* This must be included before INPUT is defined below, otherwise we
     have a name clash/problem  */
  #include <windows.h>
  #include <limits.h>
#else /* _WIN32 */
  #include <inttypes.h>
#endif /* _WIN32 */


/**
 *  CHAMELEON types
 */
typedef int8_t cham_bool_t;


/**
 * CHAMELEON Complex numbers
 */
#define CHAMELEON_HAS_COMPLEX_H 1

#if defined(_WIN32)
# include <float.h>
# if defined(__INTEL_COMPILER)
    /* Fix name conflict within the cabs prototype (_Complex) that    */
    /* conflicts with a C99 keyword.                                  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex
    #undef complex
# elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #undef  CHAMELEON_COMPLEX_CPP
    #define CHAMELEON_COMPLEX_CPP
# else
    #error "Supported compilers on WIN32 are MSVC and Intel Compiler."
# endif /* __INTEL_COMPILER */

# define isnan _isnan
# define isinf !_finite
#endif /* _WIN32 */

/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef CHAMELEON_HAS_COMPLEX_H
#endif /* __SUNPRO_C */

#ifndef __cplusplus
    #undef CHAMELEON_COMPLEX_CPP
#endif

#if defined(CHAMELEON_COMPLEX_CPP)
    #ifndef LAPACK_COMPLEX_CPP
    # define LAPACK_COMPLEX_CPP
    # warning "CHAMELEON_COMPLEX_CPP was defined, but not LAPACK_COMPLEX_CPP. Maybe you want to set both."
    #endif
    #include <complex> // needed for std::complex declaration
    #define CHAMELEON_Complex32_t std::complex<float>
    #define CHAMELEON_Complex64_t std::complex<double>
#else /* CHAMELEON_COMPLEX_CPP */
      /* not using cplusplus complex type: */

    #if defined(__STDC_NO_COMPLEX__)
    # error "Compiler support for complex number is required."
    #endif

    #define CHAMELEON_Complex32_t float  _Complex
    #define CHAMELEON_Complex64_t double _Complex

    #if CHAMELEON_HAS_COMPLEX_H
    # include <complex.h>
    #endif
#endif /* CHAMELEON_COMPLEX_CPP */

/**
 * Half precision on GPUs
 */
#if defined(__cplusplus) && defined(CHAMELEON_USE_CUDA) && (CUDA_VERSION >= 7500)
typedef __half CHAMELEON_Real16_t;
#elif defined(CHAMELEON_USE_HIP) && defined(__HIP__)
typedef hipblasHalf CHAMELEON_Real16_t;
#else
/* use short for cuda older than 7.5 and non-cuda files
 * corresponding routines would not work anyway since there is no half precision */
typedef short  CHAMELEON_Real16_t;
#endif

/**
 *  CHAMELEON Deprecated attribute
 */
#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define CHAMELEON_DEPRECATED  __attribute__((__deprecated__))
#else
#define CHAMELEON_DEPRECATED
#endif /* __GNUC__ */

BEGIN_C_DECLS

/**
 *  Global utilities
 */
static inline int chameleon_max( int a, int b ) {
    if ( a > b ) return a; else return b;
}

static inline int chameleon_min( int a, int b ) {
    if ( a < b ) return a; else return b;
}

static inline int chameleon_ceil( int a, int b ) {
    return ( a + b - 1 ) / b;
}

typedef double cham_fixdbl_t;

END_C_DECLS

#endif /* _chameleon_types_h_ */
