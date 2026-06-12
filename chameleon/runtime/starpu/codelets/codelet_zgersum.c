/**
 *
 * @file starpu/codelet_zgersum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgersum StarPU codelet
 *
 * @version 1.4.0
 * @author Romain Peressoni
 * @author Mathieu Faverge
 * @author Antoine Jego
 * @author Florent Pruvost
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgersum_redux_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( tileA->m == tileB->m );
    assert( tileA->n == tileB->n );

    TCORE_zgeadd( ChamNoTrans, tileA->m, tileA->n, 1., tileB, 1., tileA );

    (void)cl_arg;
    return;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zgersum_redux_cuda_func( void *descr[], void *cl_arg )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zone  = make_cuDoubleComplex(1.0, 0.0);
#else
    double zone  = 1.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
    cublasHandle_t  handle = starpu_cublas_get_local_handle();
    CHAM_tile_t    *tileA;
    CHAM_tile_t    *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( tileA->m == tileB->m );
    assert( tileA->n == tileB->n );

    TCUDA_zgeadd( ChamNoTrans, tileA->m, tileA->n, &zone, tileB, &zone, tileA, handle );

    (void)cl_arg;
    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgersum_redux, cl_zgersum_redux_cpu_func, cl_zgersum_redux_cuda_func, STARPU_CUDA_ASYNC )

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgersum_init_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlaset( ChamUpperLower, tileA->m, tileA->n, 0., 0., tileA );

    (void)cl_arg;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zgersum_init_cuda_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t *tileA;
    cublasStatus_t rc;

    tileA = cti_interface_get(descr[0]);

    rc = cudaMemset2D( tileA->mat, sizeof(CHAMELEON_Complex64_t) * tileA->ld, 0,
                       sizeof(CHAMELEON_Complex64_t) * tileA->m, tileA->n );
    assert( rc == CUBLAS_STATUS_SUCCESS );

    (void)cl_arg;
    (void)rc;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgersum_init, cl_zgersum_init_cpu_func, cl_zgersum_init_cuda_func, STARPU_CUDA_ASYNC )

void
RUNTIME_zgersum_set_methods( const CHAM_desc_t *A, int Am, int An )
{
#if defined(HAVE_STARPU_MPI_REDUX)
    starpu_data_set_reduction_methods( RTBLKADDR(A, ChamComplexDouble, Am, An),
                                       &cl_zgersum_redux,
                                       &cl_zgersum_init );
#endif
}

void
RUNTIME_zgersum_submit_tree( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An )
{
#if defined(HAVE_STARPU_MPI_REDUX) && defined(CHAMELEON_USE_MPI)
    starpu_mpi_redux_data_prio_tree( options->sequence->comm,
                                     RTBLKADDR(A, ChamComplexDouble, Am, An),
                                     options->priority + 1,
                                     2 /* Binary tree */ );
#else
    (void)options;
    (void)A;
    (void)Am;
    (void)An;
#endif
}

void RUNTIME_zgersum_init( void ) __attribute__( ( constructor ) );
void RUNTIME_zgersum_init( void )
{
    cl_zgersum_init.nbuffers = 1;
    cl_zgersum_init.modes[0] = STARPU_W;

    cl_zgersum_redux.nbuffers = 2;
    cl_zgersum_redux.modes[0] = STARPU_RW | STARPU_COMMUTE;
    cl_zgersum_redux.modes[1] = STARPU_R;
}
