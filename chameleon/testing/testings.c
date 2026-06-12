/**
 *
 * @file testings.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon auxiliary routines for testing structures
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @date 2025-12-19
 *
 */
#include "testings.h"

extern testing_options_t options;

/**
 * @brief List of all the testings available.
 */
static testing_t *testings = NULL;

/**
 * @brief Accuracy required to validate the tests
 */
static cham_fixdbl_t testing_accuracy = -1.;

/**
 * @brief Function to register a new testing
 */
void
testing_register( testing_t *test )
{
    assert( test->next == NULL );
    test->next = testings;
    testings = test;
}

/**
 * @brief Get the testing structure associated to a test
 */
testing_t *
testing_gettest( const char *prog_name,
                 const char *func_name )
{
    testing_t *test = testings;
    int rank = 0;
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    rank = CHAMELEON_Comm_rank();
#endif

    if ( func_name == NULL ) {
        print_usage( prog_name );
        exit(1);
    }

    while( test != NULL ) {
        /* Check the name with and without the precision */
        if ( (strcasecmp( func_name, test->name     ) == 0) ||
             (strcasecmp( func_name, test->name + 1 ) == 0) )
        {
            break;
        }
        test = test->next;
    }

    if ( test == NULL ) {
        if ( rank == 0 ) {
            printf( "Operation %s not found\n", func_name );
            printf( "The available operations are:\n" );
            test = testings;
            while( test != NULL ) {
                printf( "  %-10s %s\n", test->name, test->helper );
                test = test->next;
            }
        }
        exit(1);
    }

    return test;
}

/**
 * @brief Initialize the option structure to avoid looking for the parameters at
 * each iteration
 */
void
testing_options_init( testing_options_t *options )
{
    options->human     = parameters_getvalue_int( "human" );
    options->niter     = parameters_getvalue_int( "niter" );
    options->nowarmup  = parameters_getvalue_int( "nowarmup" );
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    options->api       = parameters_getvalue_int( "api" );
    options->async     = parameters_getvalue_int( "async" );
    options->check     = parameters_getvalue_int( "check" );
    options->forcegpu  = parameters_getvalue_int( "forcegpu" );
    options->generic   = parameters_getvalue_int( "generic" );
    options->gpus      = parameters_getvalue_int( "gpus" );
    options->mtxfmt    = parameters_getvalue_int( "mtxfmt" );
    options->P         = parameters_getvalue_int( "P" );
    options->profile   = parameters_getvalue_int( "profile" );
    options->splitsub  = parameters_getvalue_int( "splitsub" );
    options->threads   = parameters_getvalue_int( "threads" );
    options->trace     = parameters_getvalue_int( "trace" );
#else
    options->threads   = -1;
#endif

    options->file = parameters_getvalue_str( "file" );
    options->op   = parameters_getvalue_str( "op" );

    options->run_id = 0;
}

/**
 * @brief Starts the measure for the testing.
 */
void
testing_start( testdata_t *tdata )
{
#if defined(CHAMELEON_TESTINGS_VENDOR)
    /*
     * If we test the vendor functions, we want to use all the threads of the
     * machine and in order to do so we need to free the thread of the Chameleon
     * scheduler.
     */
    CHAMELEON_Pause();
#else
#if defined(CHAMELEON_USE_MPI)
    CHAMELEON_Distributed_start();
#endif
    /*
     * Create the sequence for the asynchronous calls
     */
    if ( options.async ) {
        CHAMELEON_Sequence_Create( &(tdata->sequence) );
    }

    /* Start kernel statistics */
    if ( options.profile && (options.run_id >= 0) ) {
        CHAMELEON_Enable( CHAMELEON_GENERATE_STATS );
    }

    /* Start tracing */
    if ( options.trace && (options.run_id >= 0) ) {
        CHAMELEON_Enable( CHAMELEON_GENERATE_TRACE );
    }

    /*
     * Pause the task execution if we want to time separately the task
     * submission from the task execution
     */
    if ( options.splitsub ) {
        CHAMELEON_Pause();
    }
#endif

    /* Register starting time */
    tdata->tsub  = RUNTIME_get_time();
    tdata->texec = tdata->tsub;
}

/**
 * @brief Stops the measure for the testing.
 */
void
testing_stop( testdata_t *tdata, cham_fixdbl_t flops )
{
    cham_fixdbl_t t0, t1, t2, gflops;

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    /* Submission is done, we need to start the computations */
    if ( options.async ) {
        tdata->tsub = RUNTIME_get_time();
        if ( options.splitsub ) {
            CHAMELEON_Resume();
        }
        CHAMELEON_Sequence_Wait( tdata->sequence );
        CHAMELEON_Sequence_Destroy( tdata->sequence );
    }
#if defined(CHAMELEON_USE_MPI)
    CHAMELEON_Distributed_stop();
#endif
#endif
    t2 = RUNTIME_get_time();

    t0 = tdata->texec;
    t1 = tdata->tsub;
    /*
     * texec / Submission / tsub / Execution / t
     *
     * => texec = t2 - t1
     * => tsub  = t1 - t0
     */
    tdata->tsub  = t1 - t0;
    tdata->texec = t2 - t0;

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    if ( options.splitsub ) {
        tdata->texec = t2 - t1;
    }

    /* Stop tracing */
    if ( options.trace && (options.run_id >= 0) ) {
        CHAMELEON_Disable( CHAMELEON_GENERATE_TRACE );
    }

    /* Stop kernel statistics */
    if ( options.profile && (options.run_id >= 0) ) {
        CHAMELEON_Disable( CHAMELEON_GENERATE_STATS );
    }
#endif

    gflops = flops * 1.e-9 / tdata->texec;
    run_arg_add_fixdbl( tdata->args, "time", tdata->texec );
    run_arg_add_fixdbl( tdata->args, "tsub", tdata->tsub );
    run_arg_add_fixdbl( tdata->args, "gflops", ( tdata->hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

#if defined(CHAMELEON_TESTINGS_VENDOR)
    /*
     * If we test the vendor functions, we now want to resume all the threads of
     * the machine in order to do the check and potentially the next matrix
     * initialization.
     */
    CHAMELEON_Resume();
#endif
}

/**
 * @brief Return the accuracy wanted to validate the tests
 */
cham_fixdbl_t
testing_getaccuracy()
{
#if !defined(CHAMELEON_SIMULATION)
    if ( testing_accuracy < 0. ) {
        fprintf( stderr, "WARNING: Accuracy has not been initialized yet. Please call testing_setaccuracy() prior to calling testing_getaccuracy()\n" );
        assert( 0 );
    }
#endif
    return testing_accuracy;
}

/**
 * @brief Set the global accuracy wanted to validate the tests
 */
void
testing_setaccuracy( cham_fixdbl_t accuracy )
{
    if ( accuracy <= 0. ) {
        fprintf( stderr, "WARNING: accuracy must be > 0.\n" );
        return;
    }
    testing_accuracy = accuracy;
}

