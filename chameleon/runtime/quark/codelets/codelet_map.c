/**
 *
 * @file quark/codelet_map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map Quark codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2025-01-29
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"

struct quark_map_args_s {
    cham_uplo_t          uplo;
    int                  m, n;
    cham_map_operator_t *op_fcts;
    void                *op_args;
    const CHAM_desc_t   *desc[1];
};

void CORE_map_one_quark(Quark *quark)
{
    struct quark_map_args_s *qargs = NULL;
    CHAM_tile_t             *tileA;

    quark_unpack_args_2( quark, qargs, tileA );
    qargs->op_fcts->cpufunc( qargs->op_args, qargs->uplo, qargs->m, qargs->n, 1,
                             qargs->desc[0], tileA );

    free( qargs );
}

void CORE_map_two_quark(Quark *quark)
{
    struct quark_map_args_s *qargs = NULL;
    CHAM_tile_t             *tileA;
    CHAM_tile_t             *tileB;

    quark_unpack_args_3( quark, qargs, tileA, tileB );
    qargs->op_fcts->cpufunc( qargs->op_args, qargs->uplo, qargs->m, qargs->n, 2,
                             qargs->desc[0], tileA, qargs->desc[1], tileB );

    free( qargs );
}

void CORE_map_three_quark(Quark *quark)
{
    struct quark_map_args_s *qargs = NULL;
    CHAM_tile_t             *tileA;
    CHAM_tile_t             *tileB;
    CHAM_tile_t             *tileC;

    quark_unpack_args_4( quark, qargs, tileA, tileB, tileC );
    qargs->op_fcts->cpufunc( qargs->op_args, qargs->uplo, qargs->m, qargs->n, 3,
                             qargs->desc[0], tileA, qargs->desc[1], tileB,
                             qargs->desc[2], tileC );

    free( qargs );
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    struct quark_map_args_s *qargs = NULL;
    quark_option_t          *opt   = (quark_option_t*)(options->schedopt);
    size_t                   qargs_size = 0;
    int                      i;

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    qargs_size = sizeof( struct quark_map_args_s ) + sizeof( CHAM_desc_t * ) * (ndata - 1);
    qargs = malloc( qargs_size );
    qargs->uplo    = uplo;
    qargs->m       = m;
    qargs->n       = n;
    qargs->op_fcts = op_fcts;
    qargs->op_args = op_args;
    for( i=0; i<ndata; i++ ) {
        qargs->desc[i] = data[i].desc;
    }

    switch( ndata ) {
    case 1:
        QUARK_Insert_Task(
            opt->quark, CORE_map_one_quark, (Quark_Task_Flags*)opt,
            sizeof(void*), &qargs, VALUE,
            sizeof(void*), RTBLKADDR( data[0].desc, void, m, n), cham_to_quark_access( data[0].access ),
            0 );
        break;

    case 2:
        QUARK_Insert_Task(
            opt->quark, CORE_map_two_quark, (Quark_Task_Flags*)opt,
            sizeof(void*), &qargs, VALUE,
            sizeof(void*), RTBLKADDR( data[0].desc, void, m, n), cham_to_quark_access( data[0].access ),
            sizeof(void*), RTBLKADDR( data[1].desc, void, m, n), cham_to_quark_access( data[1].access ),
            0 );
        break;

    case 3:
        QUARK_Insert_Task(
            opt->quark, CORE_map_three_quark, (Quark_Task_Flags*)opt,
            sizeof(void*), &qargs, VALUE,
            sizeof(void*), RTBLKADDR( data[0].desc, void, m, n), cham_to_quark_access( data[0].access ),
            sizeof(void*), RTBLKADDR( data[1].desc, void, m, n), cham_to_quark_access( data[1].access ),
            sizeof(void*), RTBLKADDR( data[2].desc, void, m, n), cham_to_quark_access( data[2].access ),
            0 );
        break;
    }
}
