/**
 *
 * @file starpu/runtime_mpi.h
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU mpi function implementation
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _runtime_mpi_h_
#define _runtime_mpi_h_

/**
 *  Set the tag sizes
 */
#if defined(CHAMELEON_USE_MPI)

#if !defined(HAVE_STARPU_MPI_DATA_REGISTER)
static inline void starpu_mpi_data_register( starpu_data_handle_t handle, int64_t tag, int owner )
{
    starpu_data_set_rank( handle, owner );
    starpu_data_set_tag( handle, tag );
}
#endif

#else

static inline void starpu_mpi_data_register( starpu_data_handle_t handle, int64_t tag, int owner )
{
    (void)handle;
    (void)tag;
    (void)owner;
}

#endif

#endif /* _runtime_mpi_h_ */
