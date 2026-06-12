/**
 *
 * @file cblas_wrapper.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cblas header wrapper
 *
 * @version 1.4.0
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-06-08
 *
 */
#ifndef _cblas_wrapper_h_
#define _cblas_wrapper_h_

/**
 *  CBLAS requires for scalar arguments to be passed
 *        by address rather than by value
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif

#ifndef CBLAS_H  /* cblas may have been included by another mean */
#include "coreblas/cblas.h"
#endif

/**
 * CBlas enum
 */
#ifndef CBLAS_ORDER
#define CBLAS_ORDER     enum CBLAS_ORDER
#endif
#ifndef CBLAS_TRANSPOSE
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#endif
#ifndef CBLAS_UPLO
#define CBLAS_UPLO      enum CBLAS_UPLO
#endif
#ifndef CBLAS_DIAG
#define CBLAS_DIAG      enum CBLAS_DIAG
#endif
#ifndef CBLAS_SIDE
#define CBLAS_SIDE      enum CBLAS_SIDE
#endif

#endif /* _cblas_wrapper_h_ */
