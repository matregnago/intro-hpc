/**
 *
 * @file parsec/chameleon_parsec.h
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC runtime header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_parsec_h_
#define _chameleon_parsec_h_

#include "control/common.h"

#include <parsec.h>
#include <parsec/interfaces/superscalar/insert_function.h>
#include <parsec/data_dist/matrix/matrix.h>

struct chameleon_parsec_desc_s {
    parsec_data_collection_t super;
    int                      arena_index;
    CHAM_desc_t             *desc;
    parsec_data_t          **data_map;
};

typedef struct chameleon_parsec_desc_s chameleon_parsec_desc_t;

static inline int
chameleon_parsec_get_arena_index( const CHAM_desc_t *desc ) {
    return ((chameleon_parsec_desc_t *)desc->schedopt)->arena_index;
}

static inline int
chameleon_parsec_get_arena_index_ipiv( const CHAM_ipiv_t *ipiv ) {
    assert(0);
    return -1;
}

static inline int
chameleon_parsec_get_arena_index_perm( const CHAM_ipiv_t *ipiv ) {
    assert(0);
    return -1;
}

static inline int
chameleon_parsec_get_arena_index_invp( const CHAM_ipiv_t *ipiv ) {
    assert(0);
    return -1;
}

static inline int cham_to_parsec_access( cham_access_t accessA ) {
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
#define RTBLKADDR( desc, type, m, n ) ( parsec_dtd_tile_of( (parsec_data_collection_t *) ((desc)->schedopt), \
                                                            ((parsec_data_collection_t *) (desc)->schedopt)->data_key((desc)->schedopt, m, n) ))

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _chameleon_parsec_h_ */
