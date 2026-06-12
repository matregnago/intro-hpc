/**
 *
 * @file starpu/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.4.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Gwenole Lucas
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if defined(CHAMELEON_USE_MPI) && !defined(CHAMELEON_RUNTIME_SYNC)
/**
 * This temporary patch is here to protect the DAG creation and make sure the
 * SUMMA algorithms are working fine. Please check
 * https://gitlab.inria.fr/starpu/starpu/-/merge_requests/191 for further
 * informations.
 */
static struct starpu_codelet cl_zlacpy_nowhere =
{
    .where = STARPU_NOWHERE,
    .nbuffers = 2,
    .modes = {STARPU_W, STARPU_RW},
    .model = NULL
};

static inline int
insert_task_zlacpy_nowhere( starpu_data_handle_t dst_handle,
                            starpu_data_handle_t src_handle )
{
    struct starpu_task *task = starpu_task_create();
    STARPU_ASSERT(task);
    task->name = "zlacpy_nowhere";

    task->cl = &cl_zlacpy_nowhere;

    STARPU_TASK_SET_HANDLE(task, dst_handle, 0);
    STARPU_TASK_SET_HANDLE(task, src_handle, 1);

    int ret = starpu_task_submit(task);
    STARPU_ASSERT_MSG(ret != -ENODEV, "Implementation of _starpu_data_cpy is needed for this only available architecture\n");
    STARPU_ASSERT_MSG(!ret, "Task data copy failed with code: %d\n", ret);

    return 0;
}
#endif

struct cl_zlacpy_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    int displA;
    int displB;
    int lda;
    int ldb;
};

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
static void
cl_zlacpy_rectask_func( struct starpu_task *t, void *_args )
{
    struct cl_zlacpy_args_s *clargs  = (struct cl_zlacpy_args_s *)(t->cl_arg);
    rectask_args_t          *rtargs  = (rectask_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    starpu_cham_rectask_initrequest( t, &request );

    chameleon_pzlacpy( clargs->uplo,
                       rtargs->tiles[0]->mat,
                       rtargs->tiles[1]->mat,
                       rtargs->sequence, &request );

    free( rtargs );
}
#endif /* defined(CHAMELEON_USE_RECURSIVE_TASKS) */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlacpy_starpu_func(void *descr[], void *cl_arg)
{
    static const struct starpu_data_interface_ops *interface_ops = &starpu_interface_cham_tile_ops;
    const struct starpu_data_copy_methods         *copy_methods  = interface_ops->copy_methods;
    struct cl_zlacpy_args_s                       *clargs        = (struct cl_zlacpy_args_s *)cl_arg;
    (void)clargs;

    int      workerid    = starpu_worker_get_id_check();
    unsigned memory_node = starpu_worker_get_memory_node( workerid );

    void *src_interface = descr[0];
    void *dst_interface = descr[1];

    int rc;
    (void)rc;

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    rc = copy_methods->any_to_any( src_interface, memory_node,
                                   dst_interface, memory_node, NULL );
    assert( rc == 0 );
}

static void
cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    TCORE_zlacpy( clargs->uplo, clargs->m, clargs->n, tileA, tileB );
}

static void
cl_zlacpyx_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_zlacpyx( clargs->uplo, clargs->m, clargs->n, clargs->displA,
                   tileA, clargs->lda, clargs->displB, tileB, clargs->ldb );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zlacpy_cuda_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    cublasHandle_t           handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    TCUDA_zlacpy( clargs->uplo, clargs->m, clargs->n, tileA, tileB, handle );
}

static void
cl_zlacpyx_cuda_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    cublasHandle_t           handle = starpu_cublas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCUDA_zlacpyx( clargs->uplo, clargs->m, clargs->n, clargs->displA,
                   tileA, clargs->lda, clargs->displB, tileB, clargs->ldb, handle );
}
#endif /* !defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zlacpy,  cl_zlacpy_cpu_func, cl_zlacpy_cuda_func, STARPU_CUDA_ASYNC )
CODELETS( zlacpyx, cl_zlacpyx_cpu_func, cl_zlacpyx_cuda_func, STARPU_CUDA_ASYNC )
CODELETS( zlacpy_starpu, cl_zlacpy_starpu_func, cl_zlacpy_starpu_func, STARPU_CUDA_ASYNC )

static inline void
insert_task_zlacpy_on_local_node( const RUNTIME_option_t *options,
                                  starpu_data_handle_t handleA,
                                  starpu_data_handle_t handleB )
{
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_data_cpy_priority( handleB, handleA, 0, NULL, NULL, options->priority );
#else
    starpu_data_cpy_priority( handleB, handleA, 1, NULL, NULL, options->priority );
#endif
}

#if defined(CHAMELEON_USE_MPI)
static inline void
insert_task_zlacpy_on_remote_node( const RUNTIME_option_t *options,
                                   starpu_data_handle_t handleA,
                                   starpu_data_handle_t handleB )
{
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_mpi_data_cpy_priority( handleB, handleA, options->sequence->comm, 0, NULL, NULL, options->priority );
#else
    {
        int owner = starpu_mpi_data_get_rank( handleA );
        int rank  = options->sequence->myrank;

        /* Sender side */
        if ( rank == owner ) {
            insert_task_zlacpy_nowhere( handleB, handleA );
            starpu_mpi_data_cpy_priority( handleB, handleA, options->sequence->comm, 1, NULL, NULL, options->priority );
            insert_task_zlacpy_nowhere( handleB, handleA );
        }
        /* Receiver side */
        else {
            starpu_mpi_data_cpy_priority( handleB, handleA, options->sequence->comm, 1, NULL, NULL, options->priority );
        }
    }
#endif
}
#endif

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    int                          exec       = 0;
    char                        *cl_name    = "zlacpyx";
    CHAM_tile_t                 *tileA      = A->get_blktile( A, Am, An );
    CHAM_tile_t                 *tileB      = B->get_blktile( B, Bm, Bn );
    int                          is_rectask = 0;
    rectask_args_t              *rtargs     = NULL;
    struct cl_zlacpy_args_s     *clargs     = NULL;
    enum starpu_data_access_mode accessB    = STARPU_W;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R( A, Am, An );
    CHAMELEON_ACCESS_W( B, Bm, Bn );
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

#if !defined(CHAMELEON_USE_MPI) || defined(HAVE_STARPU_MPI_DATA_CPY_PRIORITY)
    /* Insert the task */
    if ( (!options->withlacpy) &&
         (uplo == ChamUpperLower) &&
         (tileA->m == m) && (tileA->n == n) &&
         (tileB->m == m) && (tileB->n == n) &&
         (displA == 0) && (displB == 0) )
    {
#if defined(CHAMELEON_USE_MPI)
        insert_task_zlacpy_on_remote_node( options,
                                           RTBLKADDR(A, ChamComplexDouble, Am, An),
                                           RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#else
        insert_task_zlacpy_on_local_node( options,
                                          RTBLKADDR(A, ChamComplexDouble, Am, An),
                                          RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#endif
        return;
    }
#endif

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( ( tileA->format & CHAMELEON_TILE_DESC ) &&
                   ( tileB->format & CHAMELEON_TILE_DESC ) );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + sizeof(CHAM_tile_t *) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority ;
        rtargs->tiles[0] = tileA;
        rtargs->tiles[1] = tileB;
        cl_name = "zlacpy_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
        clargs->uplo   = uplo;
        clargs->m      = m;
        clargs->n      = n;
        clargs->displA = displA;
        clargs->displB = displB;
        clargs->lda    = lda;
        clargs->ldb    = ldb;
    }

    if ( uplo != ChamUpperLower ) {
        accessB = STARPU_RW;
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlacpyx,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,  RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zlacpyx ),
        STARPU_NAME, cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zlacpy )
        0 );

    (void)tileA;
    (void)tileB;
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    int                          exec       = 0;
    char                        *cl_name    = "zlacpy";
    CHAM_tile_t                 *tileA      = A->get_blktile( A, Am, An );
    CHAM_tile_t                 *tileB      = B->get_blktile( B, Bm, Bn );
    int                          is_rectask = 0;
    rectask_args_t              *rtargs     = NULL;
    struct cl_zlacpy_args_s     *clargs     = NULL;
    enum starpu_data_access_mode accessB    = STARPU_W;
    (void)rtargs;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

#if !defined(CHAMELEON_USE_MPI) || defined(HAVE_STARPU_MPI_DATA_CPY_PRIORITY)
    /* Insert the task */
    if ( (!options->withlacpy) &&
         (uplo == ChamUpperLower) &&
         (tileA->m == m) && (tileA->n == n) &&
         (tileB->m == m) && (tileB->n == n) )
    {
#if defined(CHAMELEON_USE_MPI)
        insert_task_zlacpy_on_remote_node( options,
                                           RTBLKADDR(A, ChamComplexDouble, Am, An),
                                           RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#else
        insert_task_zlacpy_on_local_node( options,
                                          RTBLKADDR(A, ChamComplexDouble, Am, An),
                                          RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#endif
        return;
    }
#endif

#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
    /* Check if this is a rectask */
    is_rectask = ( ( tileA->format & CHAMELEON_TILE_DESC ) &&
                   ( tileB->format & CHAMELEON_TILE_DESC ) );
    if ( is_rectask ) {
        rtargs = malloc( sizeof(rectask_args_t) + sizeof(CHAM_tile_t *) );
        rtargs->sequence = options->sequence;
        rtargs->parent   = options->request->parent;
        rtargs->priority = options->priority ;
        rtargs->tiles[0] = tileA;
        rtargs->tiles[1] = tileB;
        cl_name = "zlacpy_rectask";
    }
#endif

    if ( is_rectask || exec ) {
        clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
        clargs->uplo   = uplo;
        clargs->m      = m;
        clargs->n      = n;
        clargs->displA = 0;
        clargs->displB = 0;
        clargs->lda    = tileA->ld;
        clargs->ldb    = tileB->ld;
    }

    if ( uplo != ChamUpperLower ) {
        accessB = STARPU_RW;
    }

    rt_starpu_insert_task(
        &cl_zlacpy,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),

        /* Task handles */
        STARPU_R, RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,  RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( zlacpy ),
        STARPU_NAME, cl_name,

        /* Recursive task management */
        INSERT_TASK_RECTASK_PARAMS( zlacpy )
        0 );
}
