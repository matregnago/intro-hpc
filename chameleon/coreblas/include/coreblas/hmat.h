/**
 *
 * @file hmat.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU hmat-oss function declaration
 *
 * @version 1.4.0
 * @author Rocio Carratala-Saez
 * @author Mathieu Faverge
 * @author Abel Calluaud
 * @date 2024-03-27
 *
 */
#ifndef _coreblas_hmat_h_
#define _coreblas_hmat_h_

#include "coreblas.h"

#if !defined( CHAMELEON_USE_HMATOSS )
#error "This file should not be included by itself"
#endif

#include <hmat/hmat.h>

BEGIN_C_DECLS

/**
 * Functions to get linearize/unlinearize hmat into/from a buffer
 */
typedef struct hmat_buffer_comm_s {
    size_t  offset;
    char   *data;
} hmat_buffer_comm_t;

// Count the size of the structure to write
static inline void
buffer_comm_size( void * buffer, size_t n, void *user_data )
{
    size_t *size = (size_t *)user_data;
    *size += n;
    (void) buffer;
}

static inline void
buffer_comm_read(void * buffer, size_t n, void *user_data)
{
    hmat_buffer_comm_t *buffer_struct = (hmat_buffer_comm_t *) user_data;
    char *buffer_read = buffer_struct->data + buffer_struct->offset;
    memcpy(buffer, buffer_read, n);
    buffer_struct->offset += n;
}

static inline void
buffer_comm_write(void * buffer, size_t n, void *user_data)
{
    hmat_buffer_comm_t *buffer_struct = (hmat_buffer_comm_t *) user_data;
    char *buffer_write = buffer_struct->data + buffer_struct->offset;
    memcpy(buffer_write, buffer, n);
    buffer_struct->offset += n;
}

#include "coreblas/hmat_z.h"
#include "coreblas/hmat_c.h"
#include "coreblas/hmat_d.h"
#include "coreblas/hmat_s.h"

END_C_DECLS

#endif /* _coreblas_hmat_h_ */
