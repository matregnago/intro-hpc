/**
 *
 * @file chameleon.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon main header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Virouleau
 * @author Lionel Eyraud-Dubois
 * @author Alycia Lisito
 * @author Loris Lucido
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_h_
#define _chameleon_h_

/* ****************************************************************************
 * CHAMELEON types and constants
 */
#include "chameleon/config.h"
#include <stdio.h>
#include "chameleon/constants.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/struct_context.h"
#include "chameleon/descriptor_helpers.h"

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/* ****************************************************************************
 * CHAMELEON runtime common API
 */
#include "chameleon/runtime.h"

/* ****************************************************************************
 * CHAMELEON Simulation mode
 */
#include "chameleon/simulate.h"

/* ****************************************************************************
 * Include LibHQR for hierarchical trees QR/LQ factorizations
 */
#include "libhqr.h"

/* ****************************************************************************
 * CHAMELEON Tasks
 */
#include "chameleon/tasks.h"

/**
 *  @brief Structure to return information on the polar decomposition
 */
typedef struct gepdf_info_s {
    int           itQR;  /**< Number of QR iterations        */
    int           itPO;  /**< Number of Cholesky iterations  */
    cham_fixdbl_t flops; /**< Total number of flops required */
} gepdf_info_t;

#include "chameleon/chameleon_z.h"
#include "chameleon/chameleon_c.h"
#include "chameleon/chameleon_d.h"
#include "chameleon/chameleon_s.h"
#include "chameleon/chameleon_zc.h"
#include "chameleon/chameleon_ds.h"

BEGIN_C_DECLS

/* ****************************************************************************
 * CHAMELEON functionnalities
 */
int CHAMELEON_map_Tile( cham_access_t         access,
                        cham_uplo_t           uplo,
                        CHAM_desc_t          *A,
                        cham_unary_operator_t op_fct,
                        void                 *op_args ) __attribute__((deprecated("Please refer to CHAMELEON_map_Tile_Async for a more complete interface")));
int CHAMELEON_map_Tile_Async( cham_access_t         access,
                              cham_uplo_t           uplo,
                              CHAM_desc_t          *A,
                              cham_unary_operator_t op_fct,
                              void                 *op_args,
                              RUNTIME_sequence_t   *sequence,
                              RUNTIME_request_t    *request ) __attribute__((deprecated("Please refer to CHAMELEON_mapv_Tile_Async for a more complete interface")));

int CHAMELEON_mapv_Tile( cham_uplo_t          uplo,
                         int                  ndata,
                         cham_map_data_t     *data,
                         cham_map_operator_t *op_fct,
                         void                *op_args );

int CHAMELEON_mapv_Tile_Async( cham_uplo_t          uplo,
                               int                  ndata,
                               cham_map_data_t     *data,
                               cham_map_operator_t *op_fct,
                               void                *op_args,
                               RUNTIME_sequence_t  *sequence,
                               RUNTIME_request_t   *request );

/* ****************************************************************************
 * CHAMELEON Functions
 */

/* Auxiliary */
int CHAMELEON_Version           (int *ver_major, int *ver_minor, int *ver_micro);
int CHAMELEON_Initialized       (void);
int __chameleon_init            (int nworkers, int ncudas);
int __chameleon_initpar         (int nworkers, int ncudas, int nthreads_per_worker);
int __chameleon_initparcomm     (int nworkers, int ncudas, int nthreads_per_worker, MPI_Comm comm);
int __chameleon_finalize        (void);
int CHAMELEON_Pause             (void);
int CHAMELEON_Resume            (void);
int CHAMELEON_Distributed_start (void);
int CHAMELEON_Distributed_stop  (void);
int CHAMELEON_Comm_size         (void);
int CHAMELEON_Comm_rank         (void);
int CHAMELEON_Lap2Desc          ( cham_uplo_t uplo, void *A, int LDA, CHAM_desc_t *descA );
int CHAMELEON_Desc2Lap          ( cham_uplo_t uplo, CHAM_desc_t *descA, void *A, int LDA );
int CHAMELEON_Distributed_start (void);
int CHAMELEON_Distributed_stop  (void);
int CHAMELEON_Distributed_size  (int *size);
int CHAMELEON_Distributed_rank  (int *rank);
int CHAMELEON_GetThreadNbr      (void);

CHAM_context_t *CHAMELEON_GetContext();

int CHAMELEON_Lapack_to_Tile( void *A, int LDA, CHAM_desc_t *descA ) __attribute__((deprecated("Please refer to CHAMELEON_Lap2Desc() instead")));
int CHAMELEON_Tile_to_Lapack( CHAM_desc_t *descA, void *A, int LDA ) __attribute__((deprecated("Please refer to CHAMELEON_Desc2Lap() instead")));

/* Descriptor */
ssize_t CHAMELEON_Element_Size( cham_flttype_t type );

int CHAMELEON_Desc_Create_User( CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp, int mb, int nb, int bsiz,
                                int lm, int ln, int i, int j, int m, int n, int p, int q,
                                blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd,
                                blkrankof_fct_t get_rankof, void* get_rankof_arg );

int CHAMELEON_Desc_Create( CHAM_desc_t **desc, void *mat, cham_flttype_t dtyp,
                           int mb, int nb, int bsiz, int lm, int ln,
                           int i, int j, int m, int n, int p, int q );

int CHAMELEON_Desc_Create_OOC_User( CHAM_desc_t **desc, cham_flttype_t dtyp,
                                    int mb, int nb, int bsiz, int lm, int ln,
                                    int i, int j, int m, int n, int p, int q,
                                    blkrankof_fct_t get_rankof, void* get_rankof_arg );
int CHAMELEON_Desc_Create_OOC( CHAM_desc_t **desc, cham_flttype_t dtyp,
                               int mb, int nb, int bsiz, int lm, int ln,
                               int i, int j, int m, int n, int p, int q );

CHAM_desc_t *CHAMELEON_Desc_Copy( const CHAM_desc_t *descin, void *mat );
CHAM_desc_t *CHAMELEON_Desc_CopyOnZero( const CHAM_desc_t *descin, void *mat );
CHAM_desc_t *CHAMELEON_Desc_SubMatrix( CHAM_desc_t *descA, int i, int j, int m, int n );

int CHAMELEON_Desc_Destroy( CHAM_desc_t **desc );
int CHAMELEON_Desc_Acquire( const CHAM_desc_t *desc );
int CHAMELEON_Desc_Release( const CHAM_desc_t *desc );
int CHAMELEON_Desc_Flush  ( CHAM_desc_t              *desc,
                            const RUNTIME_sequence_t *sequence );

int CHAMELEON_Desc_Change_Distribution( cham_uplo_t      uplo,
                                        CHAM_desc_t     *desc,
                                        blkrankof_fct_t  new_get_rankof,
                                        void*            new_get_rankof_arg );
int CHAMELEON_Desc_Change_Distribution_Async( cham_uplo_t         uplo,
                                              CHAM_desc_t        *desc,
                                              blkrankof_fct_t     new_get_rankof,
                                              void*               new_get_rankof_arg,
                                              RUNTIME_sequence_t *sequence );

void CHAMELEON_Desc_Print( const CHAM_desc_t *desc );

/* Workspaces */
int CHAMELEON_Dealloc_Workspace (CHAM_desc_t **desc);

/* Options */
int  CHAMELEON_Enable  (int option);
int  CHAMELEON_Disable (int option);
int  CHAMELEON_Set     (int param, int  value);
int  CHAMELEON_Get     (int param, int *value);
int  CHAMELEON_Set_update_progress_callback(void (*p)(int, int)) ;

/* Sequences */
int CHAMELEON_Sequence_Create  (RUNTIME_sequence_t **sequence);
int CHAMELEON_Sequence_Destroy (RUNTIME_sequence_t *sequence);
int CHAMELEON_Sequence_Wait    (RUNTIME_sequence_t *sequence);
int CHAMELEON_Sequence_Flush   (RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int  CHAMELEON_Recursive_Desc_Create( CHAM_desc_t **descptr, void *mat,
                                      cham_flttype_t dtyp, cham_rec_t rec, int rarg,
                                      int *mb, int *nb, int lm, int ln, int m, int n, int p, int q,
                                      blkaddr_fct_t get_blkaddr, blkldd_fct_t get_blkldd,
                                      blkrankof_fct_t get_rankof, void* get_rankof_arg,
                                      const char *name );

int CHAMELEON_Ipiv_Create( CHAM_ipiv_t **ipivptr,
                           cham_side_t   side,
                           int           mb,
                           int           m,
                           int           p,
                           int           np,
                           void         *data );
int CHAMELEON_Ipiv_Destroy( CHAM_ipiv_t **ipivptr );

int  CHAMELEON_Ipiv_Flush( const CHAM_ipiv_t        *ipiv,
                           RUNTIME_sequence_t *sequence );
int  CHAMELEON_Ipiv_Gather( const CHAM_ipiv_t *ipivdesc,
                            int               *ipiv,
                            int                root );
void CHAMELEON_Ipiv_Print( const CHAM_ipiv_t *ipiv );

/* Numerical helpers */
float CHAMELEON_slamch( void );
double CHAMELEON_dlamch( void );

/**
 *
 * @ingroup Control
 *
 * @brief Initialize CHAMELEON.
 *
 ******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use.
 *
 * @param[in] gpus
 *          Number of cuda devices to use.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
#if defined(CHAMELEON_SCHED_OPENMP)

#define CHAMELEON_Init( _nworkers_, _ncudas_ )          \
    __chameleon_init( (_nworkers_), (_ncudas_) );       \
    _Pragma("omp parallel")                             \
    _Pragma("omp master")                               \
    {

#define CHAMELEON_InitPar( _nworkers_, _ncudas_, _nthreads_per_worker_ )      \
    __chameleon_initpar( (_nworkers_), (_ncudas_), (_nthreads_per_worker_) ); \
    _Pragma("omp parallel")                                                   \
    _Pragma("omp master")                                                     \
    {

#define CHAMELEON_InitParComm( _nworkers_, _ncudas_, _nthreads_per_worker_, _comm_ )        \
    __chameleon_initparcomm( (_nworkers_), (_ncudas_), (_nthreads_per_worker_), (_comm_) ); \
    _Pragma("omp parallel")                                                                 \
    _Pragma("omp master")                                                                   \
    {

#define CHAMELEON_Finalize()                    \
    }                                           \
    __chameleon_finalize();

#else

#define CHAMELEON_Init( _nworkers_, _ncudas_ )          \
    __chameleon_init( (_nworkers_), (_ncudas_) );

#define CHAMELEON_InitPar( _nworkers_, _ncudas_, _nthreads_per_worker_ ) \
    __chameleon_initpar( (_nworkers_), (_ncudas_), (_nthreads_per_worker_), MPI_COMM_WORLD );

#define CHAMELEON_InitParComm( _nworkers_, _ncudas_, _nthreads_per_worker_, _comm_ ) \
    __chameleon_initparcomm( (_nworkers_), (_ncudas_), (_nthreads_per_worker_), (_comm_) );

#define CHAMELEON_Finalize()                    \
    __chameleon_finalize();

#endif

END_C_DECLS

#endif /* _chameleon_h_ */
