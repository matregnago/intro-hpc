/**
 *
 * @file descriptor_rec.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Gwenole Lucas
 * @author Lionel Eyraud-Dubois
 * @date 2025-12-19
 *
 */
#include "control/common.h"
#include "chameleon/runtime.h"

static int
chameleon_recdesc_create( const CHAM_context_t *chamctxt,
                          const char *name, CHAM_desc_t *desc, void *mat, cham_flttype_t dtyp,
                          cham_rec_t rec, int rarg, int *mb, int *nb,
                          int lm, int ln, int m, int n, int p, int q, int i0, int j0,
                          blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd,
                          blkrankof_fct_t get_rankof, void* get_rankof_arg )
{
    CHAM_desc_t *tiledesc;
    CHAM_tile_t *tile;
    char        *subname;
    int          tempmm, tempnn;
    int          rc, i, j;

    /* Let's make sure we have at least one couple (mb, nb) defined */
    assert( (mb[0] > 0) && (nb[0] > 0) );

    /* Create the current layer descriptor */
    rc = chameleon_desc_init( chamctxt, desc, name, mat, dtyp, mb[0], nb[0],
                              lm, ln, m, n, p, q,
                              get_blkaddr, get_blkldd, get_rankof, get_rankof_arg );
    if ( rc != CHAMELEON_SUCCESS ) {
        return rc;
    }

    /* Move to the next tile size to recurse */
    mb++;
    nb++;
    if ( (mb[0] <= 0) || (nb[0] <= 0) ) {
        return CHAMELEON_SUCCESS;
    }

    for ( n=0; n<desc->nt; n++ ) {
        j = j0 * desc->nt + n; /* Used when rec = diag */

        for ( m=0; m<desc->mt; m++ ) {
            i = i0 * desc->mt + m; /* Used when rec = diag */

            switch (rec) {
            case ChamRecFull:
                break;
            case ChamRecRandom:
                /*
                 * rarg = 1: equivalent to ChamRecFull
                 * rarg = n: every n-th tiles are partitionned
                 */
                if ( random() % rarg ) continue;
                break;
            case ChamRecDiag:
                /*
                 * rarg = n: the first n tiles under and above the diagonal will
                 * be partitionned
                 */
                if ( abs( i - j ) > rarg ) continue;
                break;
            case ChamRecSmart:
                /*
                 * rarg defines the number of tiles that needs to be partitionned
                 */
                if ( n*desc->mt + m >= rarg ) continue;
                break;
            default:
                return CHAMELEON_ERR_UNEXPECTED;
            }

            tile   = desc->get_blktile( desc, m, n );
            tempmm = desc->get_blkdim( desc, m, DIM_m, desc->m );
            tempnn = desc->get_blkdim( desc, n, DIM_n, desc->n );

            chameleon_asprintf( &subname, "%s[%d,%d]", name, m, n );

            tiledesc = (CHAM_desc_t*)malloc(sizeof(CHAM_desc_t));
            rc = chameleon_recdesc_create( chamctxt, subname, tiledesc, tile->mat, desc->dtyp,
                                           rec, rarg, mb, nb,
                                           tile->ld, tempnn, /* Abuse as ln is not used */
                                           tempmm, tempnn,
                                           1, 1,             /* can recurse only on local data */
                                           i, j,             /* Used when rec = diag */
                                           chameleon_getaddr_cm, chameleon_getblkldd_cm,
                                           NULL, NULL );
            free( subname );

            tile->format = CHAMELEON_TILE_DESC;
            tile->mat    = tiledesc;

            if ( rc != CHAMELEON_SUCCESS ) {
                return rc;
            }
        }
    }

    return CHAMELEON_SUCCESS;
}

int
CHAMELEON_Recursive_Desc_Create( CHAM_desc_t **descptr, void *mat, cham_flttype_t dtyp,
                                 cham_rec_t rec, int rarg, int *mb, int *nb,
                                 int lm, int ln, int m, int n, int p, int q,
                                 blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd,
                                 blkrankof_fct_t get_rankof, void* get_rankof_arg,
                                 const char *name )
{
    CHAM_context_t *chamctxt;
    CHAM_desc_t *desc;
    int status;

    /*
     * The first layer must be allocated, otherwise we will give unitialized
     * pointers to the lower layers
     */
    assert( (mat != CHAMELEON_MAT_ALLOC_TILE) &&
            (mat != CHAMELEON_MAT_OOC) );

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_error( "CHAMELEON_Recursive_Desc_Create", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ( chamctxt->scheduler != RUNTIME_SCHED_STARPU ) {
        chameleon_error( "CHAMELEON_Recursive_Desc_Create", "CHAMELEON Recursive descriptors only available with StaRPU" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Create the current layer descriptor */
    desc = (CHAM_desc_t*)malloc(sizeof(CHAM_desc_t));
    if (desc == NULL) {
        chameleon_error("CHAMELEON_Recursive_Desc_Create", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    status = chameleon_recdesc_create( chamctxt, name, desc, mat, dtyp,
                                       rec, rarg, mb, nb,
                                       lm, ln, m, n, p, q, 0, 0,
                                       get_blkaddr, get_blkldd,
                                       get_rankof, get_rankof_arg );

    *descptr = desc;
    return status;
}
