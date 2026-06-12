/**
 *
 * @file chameleon_ztesting.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Emmanuel Agullo
 * @author Alycia Lisito
 * @author Philippe Swartvagher
 * @author Lucas Nesi
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @author Matteo Marcos
 * @author Samuel Thibault
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "testings.h"
#include "coreblas/lapacke.h"
#include <chameleon/getenv.h>

testing_options_t options;

/**
 * @brief Defines all the parameters of the testings
 *
 * @warning: keep in sync with vendor_ztesting. Must be kept in the main file for precision generation.
 * The `#if !defined(CHAMELEON_TESTINGS_VENDOR)` are not necessary in this
 * file, but they allow to easily identify which parameters are not available
 * in the vendor version.
 */
parameter_t parameters[] = {
    /* Name, helper, shname, flags, has_arg, psize, valtype, value, vallist, read, sprint */
    { "id", "Id of the run", 0, PARAM_OUTPUT, 0, 3, TestValInt, {0}, NULL, NULL, sprint_int },

    { NULL, "Options", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "help",     "Show this help",                           'h', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "human",    "Enable human readable mode",               'H', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "niter",    "Perform multiple iteration per test",      'l', PARAM_OPTION, 1, 0, TestValInt, {1}, NULL, pread_int, sprint_int },
    { "nowarmup", "Disable the warmup run to load libraries", 'w', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "check",    "Enable checking of the result",            'c', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "trace",    "Enable the trace generation",              'T', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int }, /* "--trace" and "-T" are hardcoded in a warning message below */
    { "mtxfmt",   "Change the way the matrix is stored (0: global, 1: tiles, 2: OOC)", 'F', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 1, 6, TestValInt, {1}, NULL, pread_int, sprint_int },
    { "profile",  "Display the kernel profiling",             -33, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "forcegpu", "Force kernels on GPU",                     -34, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "async",    "Switch to the Async interface",                        's', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "splitsub", "Split the task submission and execution stages",       'S', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "generic",  "Switch to the non optimized generic algorithms",       -35, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "api",      "Select the API to test (0: Descriptors, 1: Standard, 2: Lapack)", -36, PARAM_OPTION, 1, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { NULL, "Machine parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "threads", "Number of CPU workers per node",      't', PARAM_OPTION | PARAM_OUTPUT, 1, 7, TestValInt, {-1}, NULL, pread_int,    sprint_int    },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "gpus",    "Number of GPU workers per node",      'g', PARAM_OPTION | PARAM_OUTPUT, 1, 4, TestValInt, { 0}, NULL, pread_int,    sprint_int    },
    { "P",       "Rows (P) in the PxQ process grid",    'P', PARAM_OPTION | PARAM_OUTPUT, 1, 2, TestValInt, { 1}, NULL, pread_int,    sprint_int    },
    { "Q",       "Columns (Q) in the PxQ process grid", 'Q', PARAM_OUTPUT,                1, 2, TestValInt, { 1}, NULL, pread_int,    sprint_int    },
    { "custom",  "Name of custom distribution file",    -61, PARAM_OPTION,                1, 1, TestString, { 0}, NULL, pread_string, sprint_string },
#endif

    { NULL, "Main input parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "op",   "Operation to test/time ('-o help' to get the list of functions)", 'o', PARAM_OPTION | PARAM_OUTPUT, 1, 1, TestString, {0}, NULL, pread_string, sprint_string },
    { "file", "Input file",                                                      'f', PARAM_OPTION,                1, 1, TestString, {0}, NULL, pread_string, sprint_string },

    { NULL, "Matrix definition parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "m",    "Dimension M of the operation",    'm', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "n",    "Dimension N of the operation",    'n', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "k",    "Dimension K of the operation",    'k', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "nrhs", "Dimension NRHS of the operation", 'r', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "nb", "Tile size nb",       'b', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { "ib", "Inner tile size ib", 'i', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 2, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "rec", "Algorithm used to recursively tile the matrix", 'R', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 6, TestRec,    {0}, NULL, pread_rec, sprint_rec },
    { "rarg", "Parameter for the aforementioned algorithm",   '0', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "l1",  "Size of the first level of recursion",          '1', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "l2",  "Size of the second level of recursion",         '2', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "l3",  "Size of the third level of recursion",          '3', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { "lda", "Leading dimension of the matrix A", 'A', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldb", "Leading dimension of the matrix B", 'B', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldc", "Leading dimension of the matrix C", 'C', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "seedA", "Seed for the matrix A random generation", 'X', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedB", "Seed for the matrix B random generation", 'Y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedC", "Seed for the matrix C random generation", 'Z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "Matrix generation numerical parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "bump",  "Bump value to make a matrix diagonal dominant",           'z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "mode",  "Mode that specifies the eigen/singular values in xlatms", -40, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2,  4, TestValInt,       {0}, NULL, pread_int,       sprint_int       },
    { "cond",  "Conditional number of the matrix used by xlatms",         -41, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValDouble,    {0}, NULL, pread_double,    sprint_double    },

#if defined(PRECISION_z) || defined(PRECISION_d)
    { NULL, "Mixed precision Options", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "appaccuracy", "Application requested accuracy", -60, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 1, 13, TestValFixdbl, {0.}, NULL, pread_fixdbl, sprint_fixdbl },
#endif

    { NULL, "Operation specific parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "trans",  "Value of the trans parameter ('ConjTrans', 'Trans', 'NoTrans')",  -11, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transA", "Value of the transA parameter ('ConjTrans', 'Trans', 'NoTrans')", -12, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transB", "Value of the transB parameter ('ConjTrans', 'Trans', 'NoTrans')", -13, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "uplo",   "Value of the uplo parameter ('Upper', 'Lower', 'UpperLower')",    -14, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestUplo,     {0}, NULL, pread_uplo,  sprint_uplo  },
    { "diag",   "Value of the diag parameter ('NonUnit', 'Unit')",                 -15, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestDiag,     {0}, NULL, pread_diag,  sprint_diag  },
    { "side",   "Value of the side parameter ('Left', 'Right')",                   -16, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestSide,     {0}, NULL, pread_side,  sprint_side  },
    { "norm",   "Value of the norm parameter ('One', 'Frobenius', 'Inf', 'Max')",  -17, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestNormtype, {0}, NULL, pread_norm,  sprint_norm  },
    { "dir",    "Value of the dir parameter ('Forward', 'Backward')",              -18, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestDir,      {0}, NULL, pread_dir,   sprint_dir   },

    { NULL, "Operation specific scalar", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "alpha", "Value of the scalar alpha",                       'x', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "beta",  "Value of the scalar beta",                        'y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    { NULL, "QR/LQ parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "qra",    "Size of TS domain (=RH for householder trees)",           -20, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "qrp",    "Size of high level tree for distributed",                 -21, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "llvl",   "Tree used for low level reduction insides nodes",         -22, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "hlvl",   "Tree used for high level reduction between nodes",        -23, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "domino", "Enable/Disable the domino between upper and lower trees", -24, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 6, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "SVD parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "jobu",  "Value of the jobu parameter ('NoVec', 'Vec', 'Ivec', 'AllVec', 'SVec', 'OVec')",  -50, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestJob, {0}, NULL, pread_job, sprint_job },
    { "jobvt", "Value of the jobvt parameter ('NoVec', 'Vec', 'Ivec', 'AllVec', 'SVec', 'OVec')", -51, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestJob, {0}, NULL, pread_job, sprint_job },

    { NULL, "LASWP parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "k1",    "Index of the first element to permute",             -70, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "k2",    "Index of the last element to permute",              -71, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
#endif

    { "tsub",          "Graph submission time in s",             999, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "time",          "Time in s",                             1000, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "gflops",        "GFlop/s",                               1001, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "RETURN",        "Result of the testing: SUCCESS/FAILED", 1002, PARAM_OUTPUT, 2,  7, TestValInt,    {0}, NULL, pread_int,    sprint_check  },
    { "||Ax-b||",      "Norm of the residual",                  1003, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A-fact(A)||", "Norm of the residual",                  1004, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A||",         "Norm of the matrix A",                  1005, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||B||",         "Norm of the matrix B",                  1006, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||C||",         "Norm of the matrix C",                  1007, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||R||",         "Residual norm",                         1008, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||b||",         "Norm of the vector b",                  1009, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||x||",         "Norm of the vector x",                  1010, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||Ax-b||/N/eps/(||A||||x||+||b||", "",                   1011, PARAM_OUTPUT, 2, 22, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||I-QQ'||",     "Orthonormality of Q",                   1012, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },

    /* End of the list */
    { NULL, NULL, 0, 0, 0, 0, 0, {0}, NULL, NULL, NULL },
};

int main (int argc, char **argv) {

    int i;
    int rc, info = 0;
    char *input_file;
    run_list_t *runlist;
    testing_t * test;
    run_list_elt_t *run, *next;

    /* Reads the arguments from command line */
    parameters_parser( argc, argv );
    input_file = parameters_getvalue_str( "file" );
    if ( input_file != NULL ) {
        parameters_read_file( input_file );
        free(input_file);
    }

    testing_options_init( &options );

    /* Let's initialize the accuracy for the checks */
    {
#if (defined(PRECISION_z) || defined(PRECISION_d))
        cham_fixdbl_t accuracy = parameters_getvalue_fixdbl( "appaccuracy" );
        if ( accuracy > 0 ) {
            testing_setaccuracy( accuracy );
        }
        else
#endif
        {
            testing_setaccuracy( CHAMELEON_dlamch() );
        }
    }

    rc = CHAMELEON_Init( options.threads, options.gpus );
    if ( rc != CHAMELEON_SUCCESS ) {
        fprintf( stderr, "CHAMELEON_Init failed and returned %d.\n", rc );
        info = 1;
        goto end;
    }

    /* Set ncores to the right value */
    {
        parameter_t *param;
        param = parameters_get( 't' );
        param->value.ival = CHAMELEON_GetThreadNbr();

        if ( ( options.threads != -1 ) &&
             ( options.threads != param->value.ival ) )
        {
            fprintf( stderr, "WARNING: A runtime environment variable or a configuration file has overridden the number of required threads from %d to %d\n",
                     options.threads, param->value.ival );
        }
        options.threads = param->value.ival;
    }

    /* Binds the right function to be called and builds the parameters combinations */
    test = testing_gettest( argv[0], options.op );
    test_fct_t fptr = (options.api == 0) ? test->fptr_desc : test->fptr_std;
    if ( fptr == NULL ) {
        fprintf( stderr, "The %s API is not available for function %s\n",
                 (options.api == 0) ? "descriptor" : "standard", options.op );
        info = 1;
        goto end;
    }

    /* Generate the cartesian product of the parameters */
    runlist = run_list_generate( test->params );

    /* Executes the tests */
    run_print_header( test, options.check, options.human );
    run = runlist->head;

    /* Force all possible kernels on GPU */
    if ( options.forcegpu ) {
        if ( options.gpus == 0 ) {
            fprintf( stderr,
                     "--forcegpu can't be enable without GPU (-g 0).\n"
                     "  Please specify a larger number of GPU or disable this option\n" );
            info = 1;
            goto end;
        }
#if defined(CHAMELEON_SCHED_STARPU)
        int restriction = 0;
#if defined(CHAMELEON_USE_CUDA)
        restriction |= RUNTIME_CUDA;
#endif
#if defined(CHAMELEON_USE_HIP)
        restriction |= RUNTIME_HIP;
#endif
#if defined(CHAMELEON_USE_CUDA) || defined(CHAMELEON_USE_HIP)
        if ( restriction != 0 ) {
            RUNTIME_zlocality_allrestrict( restriction );
        }
#endif
        (void)restriction;
#endif
    }

    if ( !options.trace )
    {
        CHAM_context_t *chamctxt = CHAMELEON_GetContext();
        if ( (chamctxt->scheduler == RUNTIME_SCHED_STARPU) && (chameleon_getenv_get_value_int("STARPU_FXT_TRACE", 0) == 1) )
        {
            /* People often forget the --trace option */
            fprintf( stderr, "STARPU_FXT_TRACE is set, but you did not enabled Chameleon trace option with -T or --trace\n" );
        }
        /* Add here similar test for other runtime systems if necessary */
    }

    /* Warmup */
    if ( !options.nowarmup ) {
        run_arg_list_t copy = run_arg_list_copy( &(run->args) );

        /* Run the warmup test as -1 */
        options.run_id = -1;
        fptr( &copy, options.check );
        run_arg_list_destroy( &copy );
        options.run_id++;
    }

    if ( options.generic ) {
        CHAMELEON_Enable( CHAMELEON_GENERIC );
    }

    /* Perform all runs */
    while ( run != NULL ) {
        for(i=0; i<options.niter; i++) {
            run_arg_list_t copy = run_arg_list_copy( &(run->args) );
            rc = fptr( &copy, options.check );

            /* If rc < 0, we skipped the test */
            if ( rc >= 0 ) {
                run_arg_add_int( &copy, "RETURN", rc );
                run_print_line( test, &copy, options.check, options.human, options.run_id );
                options.run_id++;
                info += rc;
            }
            run_arg_list_destroy( &copy );
        }

        /* Move to next run */
        next = run->next;
        run_list_destroy( run );
        run = next;
    }

    /* Display kernel statistics if asked */
    if ( options.profile ) {
        RUNTIME_kernelprofile_display();
    }
    free( runlist );

  end:
    /* OpenMP end */
    free( options.op );
    CHAMELEON_Finalize();
    parameters_destroy();

    return info;
}
