/**
 *
 * @file run_list.c
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
 * @author Philippe Swartvagher
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "testings.h"

/**
 ********************************************************************************
 *
 * @brief Searches for a specific value by its name.
 *
 *******************************************************************************
 *
 * @param[in] arglist
 *          The list of arguments.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @retval The argument structure of the argument, NULL if not found.
 *
 *******************************************************************************
 */
const run_arg_t *
run_arg_get_byname( const run_arg_list_t *arglist, const char *name )
{
    const run_arg_t *arg = arglist->head;

    while( arg != NULL ) {
        if ( strcasecmp( name, arg->param->name ) == 0 ) {
            return arg;
        }
        arg = arg->next;
    }

    return arg;
}

/**
 ********************************************************************************
 *
 * @brief Searches for a specific value by its name.
 *
 *******************************************************************************
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 *
 *******************************************************************************
 */
val_t
run_arg_get( run_arg_list_t *arglist, const char *name, val_t defval )
{
    run_arg_t *arg = arglist->head;

    while( arg != NULL ) {
        if ( strcasecmp( name, arg->param->name ) == 0 ) {
            return arg->value;
        }
        arg = arg->next;
    }

    assert( arg  == NULL );
    arg = calloc( 1, sizeof(run_arg_t) );

    arg->param = parameters_getbyname( name );
    if ( arg->param == NULL ) {
        fprintf( stderr, "Argument %s is not registered\n", name );
        exit(1);
    }
    arg->value = defval;

    if( arglist->head == NULL ) {
        assert( arglist->tail == NULL );
        arglist->head = arg;
        arglist->tail = arg;
    }
    else {
        assert( arglist->tail != NULL );
        assert( arglist->tail->next == NULL );
        arglist->tail->next = arg;
        arglist->tail = arg;
    }

    return defval;
}

/**
 ********************************************************************************
 *
 * @brief Adds a single generic argument value by name
 *
 *******************************************************************************
 *
 * @param[inout] arglist
 *          The list of arguments to update.
 *          On exit, the argument _name_ of value _value_ is added to the list.
 *
 * @param[in] name
 *          The name of the argument to add in the list.
 *
 * @param[in] value
 *          The value of the argument to add.
 *
 * @retval 0 for success, -1 if it fails to find a paramater named _name_.
 *
 *******************************************************************************
 */
int
run_arg_add( run_arg_list_t *arglist, const char *name, val_t value )
{
    run_arg_t *arg;

    assert( arglist );

    arg = calloc( 1, sizeof(run_arg_t) );
    arg->param = parameters_getbyname( name );
    if ( arg->param == NULL ) {
        fprintf( stderr, "Argument %s does not exist\n", name );
        free( arg );
        return -1;
    }
    arg->value = value;

    if( arglist->head == NULL ) {
        assert( arglist->tail == NULL );
        arglist->head = arg;
        arglist->tail = arg;
    }
    else {
        assert( arglist->tail != NULL );
        assert( arglist->tail->next == NULL );
        arglist->tail->next = arg;
        arglist->tail = arg;
    }

    return 0;
}

/**
 * @brief Adds a single int argument value by name
 *
 * @param[inout] arglist
 *          The list of arguments to update.
 *          On exit, the argument _name_ of value _value_ is added to the list.
 *
 * @param[in] name
 *          The name of the argument to add in the list.
 *
 * @param[in] value
 *          The value of the argument to add.
 *
 * @retval 0 for success, -1 if it fails to find a paramater named _name_.
 */
int
run_arg_add_int( run_arg_list_t *arglist, const char *name, int value )
{
    val_t v;
    v.ival = value;
    return run_arg_add( arglist, name, v );
}

/**
 * @brief Adds a single double argument value by name
 *
 * @param[inout] arglist
 *          The list of arguments to update.
 *          On exit, the argument _name_ of value _value_ is added to the list.
 *
 * @param[in] name
 *          The name of the argument to add in the list.
 *
 * @param[in] value
 *          The value of the argument to add.
 *
 * @retval 0 for success, -1 if it fails to find a paramater named _name_.
 */
int
run_arg_add_double( run_arg_list_t *arglist, const char *name, double value )
{
    val_t v;
    v.dval = value;
    return run_arg_add( arglist, name, v );
}

/**
 * @brief Searches for a specific int value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
int
run_arg_get_int( run_arg_list_t *arglist, const char *name, int defval )
{
    val_t val, rval;
    val.ival = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.ival;
}

/**
 * @brief Searches for a specific float value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
float
run_arg_get_float( run_arg_list_t *arglist, const char *name, float defval )
{
    val_t val, rval;
    val.sval = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.sval;
}

/**
 * @brief Searches for a specific double value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
double
run_arg_get_double( run_arg_list_t *arglist, const char *name, double defval )
{
    val_t val, rval;
    val.dval = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.dval;
}

/**
 * @brief Searches for a specific single complex value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
CHAMELEON_Complex32_t
run_arg_get_complex32( run_arg_list_t *arglist, const char *name, CHAMELEON_Complex32_t defval )
{
    val_t val, rval;
    val.cval = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.cval;
}

/**
 * @brief Searches for a double complex value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
CHAMELEON_Complex64_t
run_arg_get_complex64( run_arg_list_t *arglist, const char *name, CHAMELEON_Complex64_t defval )
{
    val_t val, rval;
    val.zval = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.zval;
}

/**
 * @brief Searches for a cham_trans_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_trans_t
run_arg_get_trans( run_arg_list_t *arglist, const char *name, cham_trans_t defval )
{
    val_t val, rval;
    val.trans = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.trans;
}

/**
 * @brief Searches for a cham_uplo_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_uplo_t
run_arg_get_uplo( run_arg_list_t *arglist, const char *name, cham_uplo_t defval )
{
    val_t val, rval;
    val.uplo = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.uplo;
}

/**
 * @brief Searches for a cham_diag_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_diag_t
run_arg_get_diag( run_arg_list_t *arglist, const char *name, cham_diag_t defval )
{
    val_t val, rval;
    val.diag = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.diag;
}

/**
 * @brief Searches for a cham_side_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_side_t
run_arg_get_side( run_arg_list_t *arglist, const char *name, cham_side_t defval )
{
    val_t val, rval;
    val.side = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.side;
}

/**
 * @brief Searches for a cham_dir_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_dir_t
run_arg_get_dir( run_arg_list_t *arglist, const char *name, cham_dir_t defval )
{
    val_t val, rval;
    val.dir = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.dir;
}

/**
 * @brief Searches for a cham_job_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_job_t
run_arg_get_job( run_arg_list_t *arglist, const char *name, cham_job_t defval )
{
    val_t val, rval;
    val.job = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.job;
}

/**
 * @brief Searches for a cham_normtype_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_normtype_t
run_arg_get_ntype( run_arg_list_t *arglist, const char *name, cham_normtype_t defval )
{
    val_t val, rval;
    val.ntype = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.ntype;
}

/**
 * @brief Searches for the specific inner blocking value
 *
 * As opposed to other parameters, inner blocking can de defined through
 * environment variable in the Chameleon context. Thus, we first recover this
 * value as the default one, and the we get the value from the parameters.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @retval The inner blocking value.
 */
int
run_arg_get_ib( run_arg_list_t *arglist )
{
    int ib, rc;

    rc = CHAMELEON_Get( CHAMELEON_INNER_BLOCK_SIZE, &ib );
    if ( rc != CHAMELEON_SUCCESS ) {
        ib = 48;
    }
    ib = run_arg_get_int( arglist, "ib", ib );

    return ib;
}

/**
 * @brief Searches for the specific tile size value
 *
 * As opposed to other parameters, tile size can de defined through
 * environment variable in the Chameleon context. Thus, we first recover this
 * value as the default one, and the we get the value from the parameters.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @retval The inner blocking value.
 */
int
run_arg_get_nb( run_arg_list_t *arglist )
{
    int nb, rc;

    rc = CHAMELEON_Get( CHAMELEON_TILE_SIZE, &nb );
    if ( rc != CHAMELEON_SUCCESS ) {
        nb = 320;
    }
    nb = run_arg_get_int( arglist, "nb", nb );

    return nb;
}

/**
 * @brief Searches for a cham_rec_t value by its name.
 *
 * @param[inout] arglist
 *          The list of arguments.
 *          On exit, if the argument was not in the list, the default value is
 *          stored in it.
 *
 * @param[in] name
 *          The name of the argument to look for.
 *
 * @param[in] defval
 *          The default value if no argument is found with this name. This value
 *          is added to the list if not found.
 *
 * @retval The value of the argument _name_.
 */
cham_rec_t
run_arg_get_rec( run_arg_list_t *arglist, const char *name, cham_rec_t defval )
{
    val_t val, rval;
    val.rec = defval;
    rval = run_arg_get( arglist, name, val );
    return rval.rec;
}

/**
 ********************************************************************************
 *
 * @brief Frees all of memory allocated for an argument list.
 *
 *******************************************************************************
 *
 * @param[inout] args
 *          The list of arguments to free.
 *
 *******************************************************************************
 */
void run_arg_list_destroy( run_arg_list_t *arglist )
{
    run_arg_t *arg1, *arg2;

    arg1 = arglist->head;
    while( arg1 != NULL ) {
        arg2 = arg1->next;
        free( arg1 );
        arg1 = arg2;
    }

    arglist->head = NULL;
    arglist->tail = NULL;
}

/**
 ********************************************************************************
 *
 * @brief Copy a run_arg list
 *
 *******************************************************************************
 *
 * @param[in] arglist
 *          The list of running arguments to copy.
 *
 * @param[out] tailptr
 *          If tailptr is not NULL, on exit it containes the pointer to the tail
 *          of the list.
 *
 * @return The pointer to the head of the copy
 *
 *******************************************************************************
 */
run_arg_list_t
run_arg_list_copy( const run_arg_list_t *list )
{
    run_arg_list_t   copy;
    const run_arg_t *arg = list->head;
    run_arg_t *copy_curr = NULL;

    copy.head = NULL;
    copy.tail = NULL;

    /* Compute the size */
    while( arg != NULL ) {
        copy_curr = malloc( sizeof( run_arg_t ) );
        memcpy( copy_curr, arg, sizeof( run_arg_t ) );
        if ( copy.head == NULL ) {
            copy.head = copy_curr;
            copy.tail = copy_curr;
        }
        else {
            copy.tail->next = copy_curr;
            copy.tail = copy_curr;
        }
        arg = arg->next;
    }

    return copy;
}

/**
 ********************************************************************************
 *
 * @brief Add a single run argument list to the list of runs to perform
 *
 *******************************************************************************
 *
 * @param[inout] runlist
 *          The list of all the runs to perform
 *
 * @param[in] arglist
 *          The list of running arguments to copy into a new run added to the
 *          list of run.
 *
 *******************************************************************************
 */
void
run_list_add_one( run_list_t *runlist,
                  run_arg_t  *arglist )
{
    run_arg_list_t  list = { .head = arglist, .tail = NULL };
    run_list_elt_t *run;

    run = malloc( sizeof( run_list_elt_t ) );
    run->next = NULL;
    run->args = run_arg_list_copy( &list );

    if ( runlist->head == NULL ) {
        assert( runlist->tail == NULL );
        runlist->head = run;
        runlist->tail = run;
    }
    else {
        assert( runlist->tail->next == NULL );
        runlist->tail->next = run;
        runlist->tail = run;
    }
}

/**
 ********************************************************************************
 *
 * @brief Recursive function to generate the list of runs from the cartesian
 * product of the parameter values
 *
 *******************************************************************************
 *
 * @param[in] test_params
 *          The list of parameters that are considered in the test
 *
 * @param[inout] runlist
 *          The list of all the runs generated by the cartesian product.
 *
 * @param[in] arglist
 *          The current list of running arguments.
 *
 *******************************************************************************
 */
void
run_list_generate_rec( const char **test_params,
                       run_list_t  *runlist,
                       run_arg_t   *arglist )
{
    parameter_t *param = NULL;
    int          is_invalid = 1;
    run_arg_t    runarg;
    vallist_t   *vallist;

    /* End of the recursion */
    if ( *test_params == NULL ) {
        /* Add the current run_arg list to the tests */
        run_list_add_one( runlist, arglist );
        return;
    }

    /* Let's get the parameter */
    while( is_invalid && (*test_params != NULL) )
    {
        param = parameters_getbyname( *test_params );
        test_params++;

        is_invalid = ( param == NULL ) ||
            !(param->flags & PARAM_INPUT) ||
            ( param->vallist == NULL );
    }

    if ( is_invalid ) {
        /* Let's recurse one last time to register the test */
        run_list_generate_rec( test_params, runlist, arglist );
        return;
    }

    /* Let's iterate on all values */
    vallist = param->vallist;
    runarg.param = param;
    runarg.next = arglist;
    while ( vallist != NULL ) {
        runarg.value = vallist->value;
        run_list_generate_rec( test_params, runlist, &runarg );
        vallist = vallist->next;
    }

    return;
}

/**
 ********************************************************************************
 *
 * @brief Generate the list of runs from the cartesian product of the parameter
 * values.
 *
 *******************************************************************************
 *
 * @param[in] test_params
 *          The list of parameters that are considered in the test
 *
 * @return The list of all the runs generated by the cartesian product.
 *
 *******************************************************************************
 */
run_list_t *
run_list_generate( const char **params )
{
    run_list_t *runlist = calloc( 1, sizeof(run_list_t) );
    run_list_generate_rec( params, runlist, NULL );
    return runlist;
}

/**
 ********************************************************************************
 *
 * @brief Frees the run list
 *
 *******************************************************************************
 *
 * @param[inout] run
 *          The list of run to free.
 *
 *******************************************************************************
 */
void
run_list_destroy( run_list_elt_t *run )
{
    run_arg_list_destroy( &(run->args) );
    free( run );
}

/**
 * @brief The common input parameters to all tests
 */
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *common_input[]  = { "threads", NULL };
#else
const char *common_input[]  = { "threads", "gpus", "P", "Q", NULL };
#endif
/**
 * @brief The common output parameters to all tests
 */
const char *common_output[] = { "tsub", "time", "gflops", NULL };

/**
 ********************************************************************************
 *
 * @brief Print into a string the header associated to a list of parameters
 *
 *******************************************************************************
 *
 * @param[in] list
 *          The list of parameters that will be printed.
 *
 * @param[in] human
 *          Boolean to to switch between human readable and csv outputs.
 *
 * @param[in] str
 *          Pointer to the string that will store the printing
 *
 * @return The pointer to the end of the string
 *
 *******************************************************************************
 */
char *
run_print_header_partial( const char **list, int human, char *str )
{
    parameter_t *param;
    const char **pname = list;
    int rc;

    while( *pname != NULL ) {
        if ( human ) {
            param = parameters_getbyname( *pname );
            assert( param != NULL );
            switch ( param->valtype ) {
            case TestTrans:
            case TestUplo:
            case TestDiag:
            case TestSide:
            case TestJob:
            case TestNormtype:
            case TestString:
                rc = sprintf( str, " %-*s", param->psize, *pname );
                break;
            default:
                rc = sprintf( str, " %*s", param->psize, *pname );
            }
        }
        else {
            rc = sprintf( str, ";%s", *pname );
        }
        assert( rc > 0 );
        str += rc;
        pname++;
    }
    return str;
}

/**
 ********************************************************************************
 *
 * @brief Print the header associated to a test
 *
 *******************************************************************************
 *
 * @param[in] test
 *          The test that is used.
 *
 * @param[in] check
 *          Tells if check parameters should be printed or not.
 *
 * @param[in] human
 *          Boolean to to switch between human readable and csv outputs.
 *
 *******************************************************************************
 */
void
run_print_header( const testing_t *test,
                  int check, int human )
{
    int rc, rank = CHAMELEON_Comm_rank();
    char str[2048];
    char *str_ptr = str;

    if ( rank ) {
        return;
    }

    if ( human ) {
        rc = sprintf( str_ptr, "%3s %-12s",
                      "Id", "Function" );
    }
    else {
        rc = sprintf( str_ptr, "%s;%s",
                      "Id", "Function" );
    }
    str_ptr += rc;

    /* Common input */
    str_ptr = run_print_header_partial( common_input, human, str_ptr );

    /* Specific input */
    str_ptr = run_print_header_partial( test->params, human, str_ptr );

    /* Common output */
    str_ptr = run_print_header_partial( common_output, human, str_ptr );

    /* Specific output */
    str_ptr = run_print_header_partial( test->output, human, str_ptr );

    /* Specific check output */
    if ( check ) {
        run_print_header_partial( test->outchk, human, str_ptr );
    }
    fprintf( stdout, "%s\n", str );
    return;
}

/**
 ********************************************************************************
 *
 * @brief Print into a string the data associated to a list of parameters
 *
 *******************************************************************************
 *
 * @param[in] list
 *          The list of parameters that will be printed.
 *
 * @param[in] arglist
 *          The argument list in which to find the values of the parameters.
 *
 * @param[in] human
 *          Boolean to to switch between human readable and csv outputs.
 *
 * @param[in] str
 *          Pointer to the string that will store the printing
 *
 * @return The pointer to the end of the string
 *
 *******************************************************************************
 */
char *
run_print_line_partial( const char **list, const run_arg_list_t *arglist,
                        int human, char *str )
{
    parameter_t     *param;
    const char     **pname = list;
    const run_arg_t *arg;
    val_t value;

    while( *pname != NULL ) {
        arg = run_arg_get_byname( arglist, *pname );
        if ( arg == NULL ) {
            /* Should be a common parameter */
            param = parameters_getbyname( *pname );
            assert( param != NULL );

            value = param->value;
        }
        else {
            param = arg->param;
            assert( param != NULL );
            value = arg->value;
        }

        str = param->sprint( value, human, param->psize, str );

        pname++;
    }
    return str;
}

/**
 ********************************************************************************
 *
 * @brief Print the data associated to one run of a test
 *
 *******************************************************************************
 *
 * @param[in] test
 *          The test that is used.
 *
 * @param[in] arglist
 *          The argument list in which to find the values of the parameters.
 *
 * @param[in] check
 *          Tells if check parameters should be printed or not.
 *
 * @param[in] human
 *          Boolean to to switch between human readable and csv outputs.
 *
 * @param[in] id
 *          The id of the run
 *
 *******************************************************************************
 */
void
run_print_line( const testing_t *test, const run_arg_list_t *arglist,
                int check, int human, int id )
{
    int rc, rank = CHAMELEON_Comm_rank();
    char str[2048];
    char *str_ptr = str;

    if ( rank ) {
        return;
    }

    if ( human ) {
        rc = sprintf( str_ptr, "%3d %-12s",
                      id, test->name );
    }
    else {
        rc = sprintf( str_ptr, "%d;%s",
                      id, test->name );
    }
    str_ptr += rc;

    /* Common input */
    str_ptr = run_print_line_partial( common_input, arglist, human, str_ptr );

    /* Specific input */
    str_ptr = run_print_line_partial( test->params, arglist, human, str_ptr );

    /* Common output */
    str_ptr = run_print_line_partial( common_output, arglist, human, str_ptr );

    /* Specific output */
    str_ptr = run_print_line_partial( test->output, arglist, human, str_ptr );

    /* Specific check output */
    if ( check ) {
        run_print_line_partial( test->outchk, arglist, human, str_ptr );
    }

    fprintf( stdout, "%s\n", str );
    fflush(stdout);
    return;
}
