/**
 *
 * @file chameleon_tasks.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon elementary tasks main header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_tasks_h_
#define _chameleon_tasks_h_

#include "chameleon/config.h"

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include "gpucublas.h"
#endif

#if defined(CHAMELEON_USE_HIP) && !defined(CHAMELEON_SIMULATION)
#include "gpuhipblas.h"
#endif

BEGIN_C_DECLS

/**
 * @brief Kernel enum
 *
 * Those enums are Used to apply operations on specific kernels, and or for
 * tracing/profiling.
 */
typedef enum chameleon_tasktype_e {

  TASK_MAP_ONE,
  TASK_MAP_TWO,
  TASK_MAP_THREE,

  TASK_GEMM,
  TASK_HEMM,
  TASK_HER2K,
  TASK_HERK,
  TASK_SYTRF_NOPIV,
  TASK_SYMM,
  TASK_SYR2K,
  TASK_SYRK,
  TASK_TRSM,
  TASK_TRMM,

  TASK_GELQT,
  TASK_GEQRT,
  TASK_GESSM,
  TASK_GETRF,
  TASK_GETRF_INCPIV,
  TASK_GETRF_NOPIV,
  TASK_LAUUM,
  TASK_ORMLQ,
  TASK_ORMQR,
  TASK_POTRF,
  TASK_SSSSM,
  TASK_TPLQT,
  TASK_TPMLQT,
  TASK_TPMQRT,
  TASK_TPQRT,
  TASK_TRTRI,
  TASK_TSTRF,
  TASK_UNMLQ,
  TASK_UNMQR,

  TASK_GEADD,
  TASK_LASCAL,
  TASK_LACPY,
  TASK_LATRO,
  TASK_LAG2C,
  TASK_LAG2Z,
  TASK_LANGE,
  TASK_LANHE,
  TASK_LANSY,
  TASK_LASET,
  TASK_LASET2,
  TASK_PEMV,
  TASK_PLGHE,
  TASK_PLGSY,
  TASK_PLRNT,
  TASK_TILE_ZERO,

  TASK_NBKERNELS
} cham_tasktype_t;

#define TASK_TSLQT TASK_TPLQT
#define TASK_TSMLQ TASK_TPMLQT
#define TASK_TSMQR TASK_TPMQRT
#define TASK_TSQRT TASK_TPQRT
#define TASK_TTLQT TASK_TPLQT
#define TASK_TTMLQ TASK_TPMLQT
#define TASK_TTMQR TASK_TPMQRT
#define TASK_TTQRT TASK_TPQRT

typedef int (*cham_unary_operator_t)( const CHAM_desc_t *desc,
                                      cham_uplo_t uplo, int m, int n,
                                      CHAM_tile_t *data, void *op_args );

typedef int (*cham_map_cpu_fct_t)( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                                   const CHAM_desc_t *desc, CHAM_tile_t *tile, ... );

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
typedef int (*cham_map_cuda_fct_t)( cublasHandle_t handle, void *args,
                                    cham_uplo_t uplo, int m, int n, int ndata,
                                    const CHAM_desc_t *desc, CHAM_tile_t *tile, ... );
#else
typedef void *cham_map_cuda_fct_t;
#endif

#if defined(CHAMELEON_USE_HIP) && !defined(CHAMELEON_SIMULATION)
typedef int (*cham_map_hip_fct_t)( hipblasHandle_t handle, void *args,
                                   cham_uplo_t uplo, int m, int n, int ndata,
                                   const CHAM_desc_t *desc, CHAM_tile_t *tile, ... );
#else
typedef void *cham_map_hip_fct_t;
#endif

/**
 * @brief Structure to store the operator functions on any architecture
 */
typedef struct cham_map_operator_s {
    const char         *name;        /**< Name of the operator to be used in debug/tracing mode */
    cham_map_cpu_fct_t  cpufunc;     /**< Pointer to the CPU function of the operator           */
    cham_map_cuda_fct_t cudafunc;    /**< Pointer to the CUDA/cuBLAS function of the operator   */
    cham_map_hip_fct_t  hipfunc;     /**< Pointer to the HIP function of the operator           */
    int                 synchronous; /**< Specify if the kernel must be run in synchronous mode */
} cham_map_operator_t;

/**
 * @brief Structure to store the data information in the map operation
 */
typedef struct cham_map_data_s {
    cham_access_t  access; /**< Access type to the descriptor. Must be one of ChamR, ChamW, ChamRW. */
    CHAM_desc_t   *desc;   /**< Descriptor in which the data is taken to apply the map operation.   */
} cham_map_data_t;

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args );

void INSERT_TASK_gemm( const RUNTIME_option_t *options,
                       cham_trans_t transA, cham_trans_t transB,
                       int m, int n, int k, int nb,
                       double alpha, const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       double beta, const CHAM_desc_t *C, int Cm, int Cn );

void INSERT_TASK_gemmex( const RUNTIME_option_t *options,
                         cham_trans_t transA, cham_trans_t transB,
                         int m, int n, int k, int nb,
                         double alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn );

void INSERT_TASK_hgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Real16_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                  const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Real16_t beta,  const CHAM_desc_t *C, int Cm, int Cn );

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_desc_pivot_t *pivot, int k, int h, int rank );
void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2, int mt,
                               const CHAM_ipiv_t *ipivdesc, int ipivk );

#include "chameleon/tasks_z.h"
#include "chameleon/tasks_d.h"
#include "chameleon/tasks_c.h"
#include "chameleon/tasks_s.h"
#include "chameleon/tasks_zc.h"
#include "chameleon/tasks_ds.h"

END_C_DECLS

#endif /* _chameleon_tasks_h_ */
