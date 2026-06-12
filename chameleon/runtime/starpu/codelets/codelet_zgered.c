/**
 *
 * @file starpu/codelet_zgered.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgered StarPU codelet
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

void INSERT_TASK_zgered( const RUNTIME_option_t *options,
                         double threshold, int m, int n,
                         const CHAM_desc_t *A,     int Am,  int An,
                         const CHAM_desc_t *Wnorm, int Wnm, int Wnn )
{
    CHAM_tile_t          *tileA;
    double                u_low, lnorm;
    int64_t               mm, nn;
    int                   tag = -1;
    starpu_data_handle_t *handleAin;
    starpu_data_handle_t  handleAout;

    /*
     * Collect the norm of the tile on all nodes to do the the data conversion
     * if owned, and only the new data registration if not owned
     */
    {
        starpu_data_handle_t handleNorm = RTBLKADDR( Wnorm, ChamDouble, Wnm, Wnn );
        CHAM_tile_t         *tileNorm;

#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_get_data_on_all_nodes_detached( options->sequence->comm, handleNorm );
#endif
        starpu_data_acquire( handleNorm, STARPU_R );

        tileNorm = cti_handle_get( handleNorm );
        lnorm = ((double *)(tileNorm->mat))[0];

        starpu_data_release( handleNorm );
    }

    /* Get the Input handle */
    mm = Am + (A->i / A->mb);
    nn = An + (A->j / A->nb);
    handleAin = A->schedopt;
    handleAin += ((int64_t)A->lmt) * nn + mm;

    /*
     * Lets convert the tile precision based on the following criteria:
     *
     * ||A_{i,j}||_F  < u_{high} * || A ||_F / (nt * u_{low})
     * ||A_{i,j}||_F  < u_{high} * || A ||_F / nt *  1/ u_{low}
     * ||A_{i,j}||_F  < threshold / u_{low}
     */
    tileA = A->get_blktile( A, Am, An );

#if defined(CHAMELEON_USE_MPI)
    /* Backup the MPI tag */
    if ( A->myrank == tileA->rank ) {
        tag = starpu_mpi_data_get_tag( *handleAin );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

#if defined(CHAMELEON_USE_CUDA) && (CUDA_VERSION >= 7500)
#if defined(PRECISION_d)
    if ( options->withcuda ) {
        /*
         * Check for half precision
         */
        u_low = 1.e-4;
        if ( lnorm < (threshold / u_low) )
        {
#if defined(CHAMELEON_DEBUG_GERED)
            fprintf( stderr,
                     "[%2d] Convert the tile ( %d, %d ) to half precision\n",
                     A->myrank, Am, An);
#endif
            if ( A->myrank == tileA->rank )
            {
                starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexHalf );

                rt_shm_starpu_insert_task(
                    &cl_dlag2h,
                    STARPU_VALUE,            &m, sizeof(int),
                    STARPU_VALUE,            &n, sizeof(int),
                    STARPU_R,                *handleAin,
                    STARPU_W,                 handleAout,

                    /* Common task arguments */
                    INSERT_TASK_COMMON_TASK_PARAMS( dlag2h ),
                    0 );

                starpu_data_unregister_no_coherency( *handleAin );
                *handleAin     = handleAout;
                tileA->flttype = ChamComplexHalf;
                starpu_mpi_data_register( handleAout, tag, tileA->rank );
            }
            else
            {
                tileA->flttype = ChamComplexHalf;
                if ( *handleAin != NULL )
                {
                    starpu_data_unregister_no_coherency( *handleAin );
                    *handleAin = NULL;
                }
            }
            return;
        }
    }
#endif
#endif

    /*
     * Check for single precision
     */
    u_low = CHAMELEON_slamch();
    if ( lnorm < (threshold / u_low) )
    {
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert the tile ( %d, %d ) to single precision\n",
                 A->myrank, Am, An );
#endif
        if ( A->myrank == tileA->rank )
        {
            starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexFloat );

            rt_shm_starpu_insert_task(
                &cl_zlag2c,
                STARPU_VALUE,            &m, sizeof(int),
                STARPU_VALUE,            &n, sizeof(int),
                STARPU_R,                *handleAin,
                STARPU_W,                 handleAout,

                /* Common task arguments */
                INSERT_TASK_COMMON_TASK_PARAMS( zlag2c ),
                0 );

            starpu_data_unregister_no_coherency( *handleAin );
            *handleAin     = handleAout;
            tileA->flttype = ChamComplexFloat;
            starpu_mpi_data_register( *handleAin, tag, tileA->rank );
        }
        else
        {
            tileA->flttype = ChamComplexFloat;
            if ( *handleAin != NULL )
            {
                starpu_data_unregister_no_coherency( *handleAin );
                *handleAin = NULL;
            }
        }
        return;
    }
}
