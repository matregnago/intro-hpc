/**
 *
 * @file vendor_ztesting.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief MKL CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.4.0
 * @author Alycia Lisito
 * @author Lucas Nesi
 * @author Matthieu Kuhn
 * @author Philippe Swartvagher
 * @date 2025-12-19
 * @precisions normal z -> c d s
 *
 */
#include "testings.h"
#include "coreblas/lapacke.h"

testing_options_t options;

/**
 * @brief Defines all the parameters of the testings
 *
 * @warning: keep in sync with chameleon_ztesting. Must be kept in the main file for precision generation.
 */
parameter_t parameters[] = {
    /* Name, helper, shname, flags, has_arg, psize, valtype, value, vallist, read, sprint */
    { "id", "Id of the run", 0, PARAM_OUTPUT, 0, 3, TestValInt, {0}, NULL, NULL, sprint_int },

    { NULL, "Options", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "help",     "Show this help",                           'h', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "human",    "Enable human readable mode",               'H', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "niter",    "Perform multiple iteration per test",      'l', PARAM_OPTION, 1, 0, TestValInt, {1}, NULL, pread_int, sprint_int },
    { "nowarmup", "Disable the warmup run to load libraries", 'w', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "Machine parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "threads", "Number of CPU workers per node",      't', PARAM_OPTION | PARAM_OUTPUT, 1, 7, TestValInt, {-1}, NULL, pread_int, sprint_int },

    { NULL, "Main input parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "op",   "Operation to test/time ('-o help' to get the list of functions)", 'o', PARAM_OPTION | PARAM_OUTPUT, 1, 1, TestString, {0}, NULL, pread_string, sprint_string },
    { "file", "Input file",                                                      'f', PARAM_OPTION,                1, 1, TestString, {0}, NULL, pread_string, sprint_string },

    { NULL, "Matrix definition parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "m",    "Dimension M of the operation",    'm', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "n",    "Dimension N of the operation",    'n', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "k",    "Dimension K of the operation",    'k', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "nrhs", "Dimension NRHS of the operation", 'r', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "nb", "Tile size nb",       'b', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },

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

    { NULL, "Operation specific parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "trans",  "Value of the trans parameter ('ConjTrans', 'Trans', 'NoTrans')",  -11, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transA", "Value of the transA parameter ('ConjTrans', 'Trans', 'NoTrans')", -12, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transB", "Value of the transB parameter ('ConjTrans', 'Trans', 'NoTrans')", -13, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "uplo",   "Value of the uplo parameter ('Upper', 'Lower', 'UpperLower')",   -14, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestUplo,     {0}, NULL, pread_uplo,  sprint_uplo  },
    { "diag",   "Value of the diag parameter ('NonUnit', 'Unit')",   -15, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestDiag,     {0}, NULL, pread_diag,  sprint_diag  },
    { "side",   "Value of the side parameter ('Left', 'Right')",   -16, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestSide,     {0}, NULL, pread_side,  sprint_side  },
    { "norm",   "Value of the norm parameter ('One', 'Frobenius', 'Inf', 'Max')",   -17, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestNormtype, {0}, NULL, pread_norm,  sprint_norm  },

    { NULL, "Operation specific scalar", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "alpha", "Value of the scalar alpha",                       'x', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "beta",  "Value of the scalar beta",                        'y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },

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
    testing_setaccuracy( CHAMELEON_dlamch() );

    rc = CHAMELEON_Init( options.threads, 0 );
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
    test_fct_t fptr = test->fptr_std;
    if ( fptr == NULL ) {
        fprintf( stderr, "The vendor API is not available for function %s\n", options.op );
        info = 1;
        goto end;
    }

    /* Generate the cartesian product of the parameters */
    runlist = run_list_generate( test->params );

    /* Executes the tests */
    run_print_header( test, options.check, options.human );
    run = runlist->head;

    /* Warmup */
    if ( !options.nowarmup ) {
        run_arg_list_t copy = run_arg_list_copy( &(run->args) );

        /* Run the warmup test as -1 */
        options.run_id = -1;
        fptr( &copy, options.check );
        run_arg_list_destroy( &copy );
        options.run_id++;
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

    free( runlist );

  end:
    /* OpenMP end */
    free( options.op );
    CHAMELEON_Finalize();
    parameters_destroy();

    return info;
}
