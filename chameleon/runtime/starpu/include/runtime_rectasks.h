/**
 *
 * @file starpu/runtime_rectasks.h
 *
 * @copyright 2019-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU wont use implementations to flush pieces of data
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _runtime_rectasks_h_
#define _runtime_rectasks_h_

/**
 * @brief Recursive task arguments structure
 */
struct rectask_args_s {
    /**
     * @brief Defines the context in which the tasks are submitted, it
     * will be destroyed only after a synchronization that guaranties the
     * execution of all tasks. Thus, we can hold here a pointer to sumbit rectask
     * tasks within the same context.
     */
    RUNTIME_sequence_t *sequence;
    /**
     * @brief Parent task of all the tasks that will be submited within
     * this rectask. Used for profiling information.
     */
    struct starpu_task *parent;
    /**
     * @brief Priority of the parent task to propagate the information to the
     * subtasks if needed.
     */
    int priority;
    /**
     * @brief List of tiles used to submit the tasks. Each codelet uses a
     * different amount of tiles, thus we store at the end of the structure a
     * variadic size field to store the list of the tiles.
     */
    CHAM_tile_t *tiles[1];
};
typedef struct rectask_args_s rectask_args_t;

/**
 * @brief Macro to specify the parent tasks to improve profiling informations
 */
#if defined(CHAMELEON_RECURSIVE_TASKS_PROFILE)
#define INSERT_TASK_RECTASK_PROFILE_PARAM STARPU_RECURSIVE_TASK_PARENT, (rtargs ? rtargs->parent : NULL ),
#else
#define INSERT_TASK_RECTASK_PROFILE_PARAM
#endif

/**
 * @brief Macro to insert and factorize the recursive task arguments passed to the rt_starpu_insert_task function.
 */
#if defined(CHAMELEON_USE_RECURSIVE_TASKS)
#define INSERT_TASK_RECTASK_PARAMS(__name__)                              \
    STARPU_RECURSIVE_TASK_FUNC, (is_rectask ? starpu_cham_is_recursive : starpu_cham_is_not_recursive ), \
    STARPU_RECURSIVE_TASK_GEN_DAG_FUNC,     cl_##__name__##_rectask_func, \
    STARPU_RECURSIVE_TASK_GEN_DAG_FUNC_ARG, rtargs,                       \
    INSERT_TASK_RECTASK_PROFILE_PARAM
#else
#define INSERT_TASK_RECTASK_PARAMS(__name__)
#endif

/**
 * @brief Internal function to return true if the task will always be recursive
 */
static inline int
starpu_cham_is_recursive( __attribute__((unused)) struct starpu_task *task,
                          __attribute__((unused)) void               *args,
                          __attribute__((unused)) void              **descs )
{
    return 1;
}

/**
 * @brief Internal function to return false if the task will never be recursive
 */
static inline int
starpu_cham_is_not_recursive( __attribute__((unused)) struct starpu_task *task,
                              __attribute__((unused)) void               *args,
                              __attribute__((unused)) void              **descs )
{
    return 0;
}

/**
 * @brief Internal function to initialize request in recursive tasks
 */
static inline void
starpu_cham_rectask_initrequest( struct starpu_task *task,
                                 RUNTIME_request_t  *request )
{
    /* Register the task parent */
    request->parent = task;

    /* Propagate the rectask priority to subtasks if needed */
    request->priority = task->priority;

    return;
}

#endif /* _runtime_rectasks_h_ */
