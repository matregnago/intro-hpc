/**
 *
 * @file parsec/codelet_map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map PaRSEC codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks.h"

struct parsec_map_args_s {
    cham_uplo_t          uplo;
    int                  m, n;
    cham_map_operator_t *op_fcts;
    void                *op_args;
    const CHAM_desc_t   *desc[1];
};

static inline int
CORE_map_one_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    const CHAM_desc_t        *descA;
    CHAM_tile_t               tileA;

    parsec_dtd_unpack_args( this_task, &pargs, &(tileA.mat) );

    descA = pargs->desc[0];
    tileA.rank    = 0;
    tileA.m       = (pargs->m == (descA->mt-1)) ? (descA->m - pargs->m * descA->mb) : descA->mb;
    tileA.n       = (pargs->n == (descA->nt-1)) ? (descA->n - pargs->n * descA->nb) : descA->nb;
    tileA.ld      = descA->get_blkldd( descA, pargs->m );
    tileA.format  = CHAMELEON_TILE_FULLRANK;
    tileA.flttype = descA->dtyp;

    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 1,
                             descA, &tileA );

    free( pargs );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

static inline int
CORE_map_two_parsec( parsec_execution_stream_t *context,
                     parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    const CHAM_desc_t        *descA, *descB;
    CHAM_tile_t               tileA,  tileB;

    parsec_dtd_unpack_args( this_task, &pargs, &(tileA.mat), &(tileB.mat) );

    descA = pargs->desc[0];
    tileA.rank    = 0;
    tileA.m       = (pargs->m == (descA->mt-1)) ? (descA->m - pargs->m * descA->mb) : descA->mb;
    tileA.n       = (pargs->n == (descA->nt-1)) ? (descA->n - pargs->n * descA->nb) : descA->nb;
    tileA.ld      = descA->get_blkldd( descA, pargs->m );
    tileA.format  = CHAMELEON_TILE_FULLRANK;
    tileA.flttype = descA->dtyp;

    descB = pargs->desc[1];
    tileB.rank    = 0;
    tileB.m       = (pargs->m == (descB->mt-1)) ? (descB->m - pargs->m * descB->mb) : descB->mb;
    tileB.n       = (pargs->n == (descB->nt-1)) ? (descB->n - pargs->n * descB->nb) : descB->nb;
    tileB.ld      = descB->get_blkldd( descB, pargs->m );
    tileB.format  = CHAMELEON_TILE_FULLRANK;
    tileB.flttype = descB->dtyp;

    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 2,
                             descA, &tileA, descB, &tileB );

    free( pargs );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

static inline int
CORE_map_three_parsec( parsec_execution_stream_t *context,
                       parsec_task_t             *this_task )
{
    struct parsec_map_args_s *pargs = NULL;
    const CHAM_desc_t        *descA, *descB, *descC;
    CHAM_tile_t               tileA,  tileB,  tileC;

    parsec_dtd_unpack_args( this_task, &pargs, &(tileA.mat), &(tileB.mat), &(tileC.mat) );

    descA = pargs->desc[0];
    tileA.rank    = 0;
    tileA.m       = (pargs->m == (descA->mt-1)) ? (descA->m - pargs->m * descA->mb) : descA->mb;
    tileA.n       = (pargs->n == (descA->nt-1)) ? (descA->n - pargs->n * descA->nb) : descA->nb;
    tileA.ld      = descA->get_blkldd( descA, pargs->m );
    tileA.format  = CHAMELEON_TILE_FULLRANK;
    tileA.flttype = descA->dtyp;

    descB = pargs->desc[1];
    tileB.rank    = 0;
    tileB.m       = (pargs->m == (descB->mt-1)) ? (descB->m - pargs->m * descB->mb) : descB->mb;
    tileB.n       = (pargs->n == (descB->nt-1)) ? (descB->n - pargs->n * descB->nb) : descB->nb;
    tileB.ld      = descB->get_blkldd( descB, pargs->m );
    tileB.format  = CHAMELEON_TILE_FULLRANK;
    tileB.flttype = descB->dtyp;

    descC = pargs->desc[2];
    tileC.rank    = 0;
    tileC.m       = (pargs->m == (descC->mt-1)) ? (descC->m - pargs->m * descC->mb) : descC->mb;
    tileC.n       = (pargs->n == (descC->nt-1)) ? (descC->n - pargs->n * descC->nb) : descC->nb;
    tileC.ld      = descC->get_blkldd( descC, pargs->m );
    tileC.format  = CHAMELEON_TILE_FULLRANK;
    tileC.flttype = descC->dtyp;

    pargs->op_fcts->cpufunc( pargs->op_args, pargs->uplo, pargs->m, pargs->n, 3,
                             descA, &tileA, descB, &tileB, descC, &tileC );

    free( pargs );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    parsec_taskpool_t        *PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    struct parsec_map_args_s *pargs      = NULL;
    size_t                    pargs_size = 0;
    int                       i;

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    pargs_size = sizeof( struct parsec_map_args_s ) + sizeof( CHAM_desc_t * ) * (ndata - 1);
    pargs = malloc( pargs_size );
    pargs->uplo    = uplo;
    pargs->m       = m;
    pargs->n       = n;
    pargs->op_fcts = op_fcts;
    pargs->op_args = op_args;
    for( i=0; i<ndata; i++ ) {
        pargs->desc[i] = data[i].desc;
    }

    switch( ndata ) {
    case 1:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_one_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PARSEC_DTD_ARG_END );
        break;

    case 2:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_two_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PASSED_BY_REF, RTBLKADDR( data[1].desc, void, m, n ), chameleon_parsec_get_arena_index( data[1].desc ) | cham_to_parsec_access( data[1].access ),
            PARSEC_DTD_ARG_END );
        break;

    case 3:
        parsec_dtd_taskpool_insert_task(
            PARSEC_dtd_taskpool, CORE_map_three_parsec, options->priority, op_fcts->name,
            sizeof(struct parsec_map_args_s*), &pargs, VALUE,
            PASSED_BY_REF, RTBLKADDR( data[0].desc, void, m, n ), chameleon_parsec_get_arena_index( data[0].desc ) | cham_to_parsec_access( data[0].access ),
            PASSED_BY_REF, RTBLKADDR( data[1].desc, void, m, n ), chameleon_parsec_get_arena_index( data[1].desc ) | cham_to_parsec_access( data[1].access ),
            PASSED_BY_REF, RTBLKADDR( data[2].desc, void, m, n ), chameleon_parsec_get_arena_index( data[2].desc ) | cham_to_parsec_access( data[2].access ),
            PARSEC_DTD_ARG_END );
        break;
    }
}
