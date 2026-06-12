/**
 *
 * @file starpu/codelet_zgerst.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgerst StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Ana Hourcau
 * @date 2024-10-18
 * @precisions normal z -> d
 *
 */
#include "chameleon_starpu_internal.h"
#include <coreblas/lapacke.h>
#include "runtime_codelet_zc.h"
#include "runtime_codelet_z.h"

void INSERT_TASK_zgerst( const RUNTIME_option_t *options,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t          *tileA;
    int64_t               mm, nn;
    int                   tag = -1;
    starpu_data_handle_t *handleAin;
    starpu_data_handle_t  handleAout;

    tileA = A->get_blktile( A, Am, An );

    /* Get the Input handle */
    mm = Am + (A->i / A->mb);
    nn = An + (A->j / A->nb);
    handleAin = A->schedopt;
    handleAin += ((int64_t)A->lmt) * nn + mm;

    if ( tileA->flttype == ChamComplexDouble )
    {
        starpu_data_handle_t *copy = handleAin;

        /* Remove first copy */
        copy += ((int64_t)A->lmt * (int64_t)A->lnt);
        if ( *copy ) {
            starpu_data_unregister_no_coherency( *copy );
            *copy = NULL;
        }

        /* Remove second copy */
        copy += ((int64_t)A->lmt * (int64_t)A->lnt);
        if ( *copy ) {
            starpu_data_unregister_no_coherency( *copy );
            *copy = NULL;
        }

        return;
    }

    if ( A->myrank != tileA->rank )
    {
        tileA->flttype = ChamComplexDouble;
        if ( *handleAin != NULL )
        {
            starpu_data_unregister_no_coherency( *handleAin );
            *handleAin = NULL;
        }
        return;
    }

#if defined(CHAMELEON_USE_MPI)
    tag = starpu_mpi_data_get_tag( *handleAin );
#endif /* defined(CHAMELEON_USE_MPI) */

    starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexDouble );

    switch( tileA->flttype ) {
#if defined(CHAMELEON_USE_CUDA) && (CUDA_VERSION >= 7500)
#if defined(PRECISION_d)
        /*
         * Restore from half precision
         */
    case ChamComplexHalf:
        assert( options->withcuda );
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert back the tile ( %d, %d ) from half precision\n",
                 A->myrank, Am, An );
#endif
        rt_shm_starpu_insert_task(
            &cl_hlag2d,
            STARPU_VALUE,            &m,                 sizeof(int),
            STARPU_VALUE,            &n,                 sizeof(int),
            STARPU_R,                *handleAin,
            STARPU_W,                 handleAout,

            /* Common task arguments */
            INSERT_TASK_COMMON_TASK_PARAMS( hlag2d ),
            0 );
        break;
#endif
#endif

    case ChamComplexFloat:
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert back the tile ( %d, %d ) from single precision\n",
                 A->myrank, Am, An );
#endif

        rt_shm_starpu_insert_task(
            &cl_clag2z,
            STARPU_VALUE,            &m,                 sizeof(int),
            STARPU_VALUE,            &n,                 sizeof(int),
            STARPU_R,                *handleAin,
            STARPU_W,                 handleAout,

            /* Common task arguments */
            INSERT_TASK_COMMON_TASK_PARAMS( clag2z ),
            0 );
        break;

    default:
        fprintf( stderr, "ERROR: Unknonw input datatype" );
    }

    starpu_data_unregister_no_coherency( *handleAin );
    *handleAin     = handleAout;
    tileA->flttype = ChamComplexDouble;
    starpu_mpi_data_register( handleAout, tag, tileA->rank );
}
