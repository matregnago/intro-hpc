/**
 *
 * @file parameters.c
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
 * @author Lionel Eyraud-Dubois
 * @author Lucas Barros De Assis
 * @date 2025-12-19
 *
 */
#include "testings.h"

extern parameter_t parameters[];

void print_usage( const char* prog_name )
{
    parameter_t *param = parameters;

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    if ( CHAMELEON_Initialized() ) {
        int rank = CHAMELEON_Comm_rank();
        if (rank != 0) {
            return;
        }
    }
#endif

    printf( "Usage:\n"
            "  %s -o|--op operation_name [options]\n"
            "  %s -f|--file input_file [options]\n",
            prog_name, prog_name );

    while( param->helper != NULL )
    {
        char str[STR_MAX_LENGTH];

        /* This is not an option, we skip it */
        if ( !(param->flags & PARAM_OPTION) ) {
            param++;
            continue;
        }

        /* This is an option header */
        if ( param->name == NULL ) {
            printf( "\n  %s:\n", param->helper );
            param++;
            continue;
        }

        if ( param->shname > 0 ) {
            snprintf( str, STR_MAX_LENGTH, "-%c, --%s",
                      param->shname, param->name );
        }
        else {
            snprintf( str, STR_MAX_LENGTH, "    --%s",
                      param->name );
        }

        /* If an argument is needed, add " x" */
        if ( param->has_arg > 0 ) {
            int len = strlen(str);
            assert( len < (STR_MAX_LENGTH-3) );

            str[ len   ] = ' ';
            str[ len+1 ] = 'x';
            str[ len+2 ] = '\0';
        }
        /* Parameter for all tests which has a defautl value */
        if ( param->has_arg == 1 ) {
            char defaultstr[STR_MAX_LENGTH];
            param->sprint( param->value, 1, 1, defaultstr );
            printf( "    %-23s %s [default:%s]\n",
                    str, param->helper, defaultstr );
        }
        else {
            printf( "    %-23s %s\n",
                    str, param->helper );
        }
        param++;
    }

    printf( "\n"
            "For example: %s -H -o gemm -t 2 -m 2000 -n 2000 -k 2000 -b 200\n"
            "  will run one gemm with three matrices of size 2000x2000 each and a tile size of 200.\n"
            "  The output will be in the human readable format\n"
            "\n", prog_name );
#if !defined(CHAMELEON_TESTINGS_VENDOR)
    printf( "Remarks about timing:\n"
            "  Timings are reported respectively as 'tsub' for the graph submission time, and 'time'\n"
            "  for the execution time.\n"
            "  By default the synchronous tile interface is used to perform the timings. 'tsub' is null.\n"
            "  If the --async option is enabled, then the asynchronous interface is called. 'tsub' reports\n"
            "  the task submission time, and 'time' the execution time that includes 'tsub'.\n"
            "  If the --splitsub option is enabled, then the asynchronous interface is called and task\n"
            "  submission is fully performed before starting the computation. 'tsub' reports the\n"
            "  task submission time, and 'time' the execution time excluding 'tsub'.\n"
            "  Note that the 'gflops' field is always computed with 'time'\n" );
#endif
}

/**
 ********************************************************************************
 *
 * @brief Get the list of values associated to a given parameter
 *
 *******************************************************************************
 *
 * @param[in] name
 *          The name of the parameter we are interested in.
 *
 * @return NULL if no parameter exists with this name, otherwise the pointer to
 * the list of values associated to this parameter.
 *
 *******************************************************************************
 */
vallist_t *
parameters_getlist( const char *name )
{
    parameter_t *param = parameters_getbyname( name );
    if ( param == NULL ) {
        return NULL;
    }
    else {
        return param->vallist;
    }
}

/**
 ********************************************************************************
 *
 * @brief Parses a list in form A1, A2, ..., An and insert the values in an
 * argument list.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the list of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] liststr
 *          The string that holds the list
 *
 *******************************************************************************
 */
void
parameters_read_list( parameter_t *param,
                      const char  *liststr )
{
    const char *delim = ", \n";
    char *str = strdup( liststr );
    char *token, *saveptr;
    vallist_t *previous, *current;

    /* Initialize the list items */
    previous = NULL;
    current  = param->vallist;

    /* Move to the end of the list if some parameters have already been registered */
    while( current != NULL ) {
        previous = current;
        current  = current->next;
    }
    assert( current == NULL );

    token = strtok_r( str, delim, &saveptr );
    while ( token != NULL ) {
        current = calloc( 1, sizeof(vallist_t) );

        /* Read the value */
        current->value = param->read( token );

        /* Insert at the end of the list */
        if ( previous != NULL ) {
            previous->next = current;
        }
        else {
            /* Nothing was in the list */
            param->vallist = current;
        }

        previous = current;

        /* Move to the next token */
        token = strtok_r( NULL, delim, &saveptr );
    }

    free( str );
}

/**
 ********************************************************************************
 *
 * @brief Parses a list in form start:end[:step] and inserts the values in an
 * argument list.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the range of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] rangestr
 *          The string that holds the range
 *
 * @param[in] min
 *          The minimum value available
 *
 * @param[in] max
 *          The maximum value available
 *
 *******************************************************************************
 */
void
parameters_read_intrange( parameter_t *param,
                          const char  *rangestr,
                          int min, int max )
{
    int start, end, step, count;
    vallist_t *previous, *current;

    max = (max == -1) ? INT32_MAX : max;

    count = sscanf( rangestr, "%d:%d:%d", &start, &end, &step );
    if ( count < 2 ) {
        fprintf(stderr, "Incorrect range syntax (%s): data skipped\n", rangestr );
        return;
    }
    else if (count == 2) {
        step = 1;
    }

    /* Check the range */
    if ( (start < min) || (start > max) ||
         (end   < min) || (end   > max) )
    {
        /* Try to shift to 0 to see if now we fit */
        start += min;
        end   += min;
        if ( (start < min) || (start > max) ||
             (end   < min) || (end   > max) )
        {
            fprintf( stderr, "Incorrect range values outside the possible ranges [%d:%d]",
                     min, max );
            if ( min > 0 ) {
                fprintf( stderr, " or [%d:%d]\n", 0, max-min );
            }
            else {
                fprintf( stderr, "\n" );
            }
        }
    }

    /* Initialize the list items */
    previous = NULL;
    current  = param->vallist;

    /* Move to the end of the list if some parameters have already been registered */
    while( current != NULL ) {
        previous = current;
        current  = current->next;
    }
    assert( current == NULL );

    while ( start <= end ) {
        current = calloc( 1, sizeof(vallist_t) );

        /* Read the value */
        current->value.ival = start;

        /* Insert at the end of the list */
        if ( previous != NULL ) {
            previous->next = current;
        }
        else {
            /* Nothing was in the list */
            param->vallist = current;
        }

        previous = current;

        start += step;
    }
}

/**
 ********************************************************************************
 *
 * @brief Wrapper to parse a list or range of values associated to a parameter.
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter associated to the list.
 *          On exit, the range of values are added to the parameter list of
 *          possible values.
 *
 * @param[in] values
 *          The string that holds the range of list of values
 *
 *******************************************************************************
 */
void
parameters_read( parameter_t *param,
                 const char  *values )
{
    int range;

    if ( values == NULL ) {
        fprintf( stderr, "Error passing NULL parameter to read\n" );
        exit(EXIT_FAILURE);
    }
    range = (values != NULL) && ( strchr( values, ':' ) != NULL );

    /* If we have a ranged of integer values */
    if ( range )
    {
        switch ( param->valtype ) {
        case TestValInt:
            parameters_read_intrange( param, values, 0, -1 );
            break;
        case TestTrans:
            parameters_read_intrange( param, values, ChamNoTrans, ChamConjTrans );
            break;
        case TestUplo:
            parameters_read_intrange( param, values, ChamUpper, ChamUpperLower );
            break;
        case TestDiag:
            parameters_read_intrange( param, values, ChamNonUnit, ChamUnit );
            break;
        case TestSide:
            parameters_read_intrange( param, values, ChamLeft, ChamRight );
            break;
        case TestNormtype:
            parameters_read_intrange( param, values, ChamOneNorm, ChamMaxNorm );
            break;
        default:
            fprintf( stderr, "parameters_read: range is not available for this datatype (%d)\n",
                     param->valtype );
        }
        return;
    }

    parameters_read_list( param, values );
}

/**
 ********************************************************************************
 *
 * @brief Generic function to add value(s) to a given parameter
 *
 *******************************************************************************
 *
 * @param[inout] param
 *          The parameter that will receive the value
 *          On exit, the value(s) (switch, list, range, ...) is/are added to the
 *          parameter list of possible values
 *
 * @param[in] values
 *          The string that holds the values (list, range, NULL if switch)
 *
 *******************************************************************************
 */
void
parameters_addvalues( parameter_t *param,
                      const char  *values )
{
    if ( param->has_arg == 0 ) {
        param->value.ival = 1;
    }
    else if ( param->has_arg == 1 ) {
        param->value = param->read( values );
    }
    else {
        parameters_read( param, values );
    }
}

/**
 ********************************************************************************
 *
 * @brief Parses an input test file.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The name of the input file.
 *
 *******************************************************************************
 */
void
parameters_read_file( const char  *filename )
{
    FILE        *fp;
    const char  *delim = " =";
    char        *saveptr;
    char        *line_read, *line;
    char        *name, *values;
    size_t       len = 256;
    parameter_t *param;

    fp = fopen( filename, "r" );
    if ( fp == NULL ) {
        fprintf( stderr, "Error reading input file %s\n", filename );
        perror("fopen");
        exit(1);
    }

    len = 256;
    line_read = malloc( sizeof( char ) * len );

    while ( getline( &line_read, &len, fp ) != -1 )
    {
        line = line_read;

        /* Ignores comments and empty lines */
        if ( (line[0] == '#' ) ||
             (line[0] == '\n') )
        {
            continue;
        }

        /* Removes possible extra spaces */
        while ( line[0] == ' ' ) {
            line++;
        }

        /* Reads the parameter name and values */
        name   = strtok_r( line, delim, &saveptr );
        values = strtok_r( NULL, "",    &saveptr );

        /* Goes for the listed values */
        while ( (values[0] == ' ') ||
                (values[0] == '=') )
        {
            values++;
        }

        param = parameters_getbyname( name );
        if ( param == NULL ) {
            fprintf( stderr, "Parameter %s is not know. We skip it\n", name );
            continue;
        }
        parameters_addvalues( param, values );
    }

    free(line_read);
    fclose(fp);
}

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
parameters_compute_q( int p )
{
    parameter_t *param;
    int np = CHAMELEON_Comm_size();

    if ( (np % p) != 0 ) {
        fprintf( stderr, "ERROR: The number of processes (%d) must be a multiple of P (%d)\n", np, p );
        exit(1);
    }

    param = parameters_get( 'Q' );
    param->value.ival = np / p;
    return param->value.ival;
}

/**
 ********************************************************************************
 *
 * @brief Helper function to generate the testing descriptors with the right
 * data distrbution.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The name of the input file.
 *
 *******************************************************************************
 */
int
parameters_desc_create( const char *id, CHAM_desc_t **descptr, cham_flttype_t dtyp,
                        int mb, int nb, int lm, int ln, int m, int n )
{
    custom_dist_t *custom_args = NULL;
    const char    *custom = parameters_getvalue_str( "custom" );
    intptr_t       mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int            rc;

    mtxfmt = -mtxfmt; /* Inverse sign to get the defined values */

    if ( cham_is_mixed( dtyp )  && ( (void*)mtxfmt != CHAMELEON_MAT_ALLOC_TILE ) )
    {
        fprintf( stderr, "parameters_desc_create: Mixed precision descriptors can only be used with tiled allocation. Please enforce '--mtxfmt=1' \n" );
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    if ( !custom ) {
        int P = parameters_getvalue_int( "P" );
        int Q = parameters_compute_q( P );
        rc = CHAMELEON_Desc_Create(
            descptr, (void*)mtxfmt, dtyp, mb, nb, mb * nb, lm, ln, 0, 0, m, n, P, Q );
        /* Dirty hack to name the decriptor as wanted by the testing, should be modified in the future */
        if ( (*descptr)->name ) {
            free( (*descptr)->name );
        }
        (*descptr)->name = strdup(id);
        return rc;
    }

    if ( ((void*)mtxfmt) == CHAMELEON_MAT_ALLOC_GLOBAL ) {
        fprintf( stderr, "In parameters_desc_create, cannot use custom distributions with global matrix allocation (Use --mtxfmt=1)\n" );
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    rc = chameleon_getrankof_custom_init( &custom_args, custom );
    if ( rc != CHAMELEON_SUCCESS ) {
        return rc;
    }

    rc = CHAMELEON_Desc_Create_User(
        descptr, (void*)mtxfmt, dtyp, mb, nb, mb * nb, lm, ln, 0, 0, m, n, CHAMELEON_Comm_size(), 1,
        NULL, NULL, chameleon_getrankof_custom, custom_args );
    /* Dirty hack to name the decriptor as wanted by the testing, should be modified in the future */
    if ( (*descptr)->name ) {
        free( (*descptr)->name );
    }
    (*descptr)->name = strdup(id);
    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Helper function to destroy the testing descriptors.
 *
 *******************************************************************************
 *
 * @param[inout] descptr
 *      The descriptor to destroy. On exit the descriptor can no longer be used.
 *
 *******************************************************************************
 */
int
parameters_desc_destroy(CHAM_desc_t **descptr)
{
    CHAM_desc_t *desc;
    if ( descptr == NULL ) {
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    desc = *descptr;
    if ( desc == NULL ) {
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( desc->get_rankof_init_arg ) {
        if ( desc->get_rankof_init == chameleon_getrankof_custom ) {
            chameleon_getrankof_custom_destroy( (custom_dist_t**)&(desc->get_rankof_init_arg) );
        }
    }
    return CHAMELEON_Desc_Destroy( descptr );
}
#endif

void
parameters_getopt_init( char           *optstring,
                        struct option **longopts )
{
    parameter_t *param = parameters;
    int i;
    int nboptions = 0;
    int strpos = 0;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( !(param->flags & PARAM_OPTION) ||
             (param->name == NULL) )
        {
            param++;
            continue;
        }

        nboptions++;

        if ( param->shname < 0 ) {
            param++;
            continue;
        }

        optstring[strpos] = param->shname;
        strpos++;
        assert( strpos < STR_MAX_LENGTH );

        if ( param->has_arg > 0 ) {
            optstring[strpos] = ':';
            strpos++;
            assert( strpos < STR_MAX_LENGTH );
        }
        param++;
    }
    optstring[strpos] = '\0';

    /* Now, let's generate the long opt if needed */
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    if ( longopts != NULL ) {
        struct option *opt;
        *longopts = calloc( nboptions+1, sizeof( struct option ) );

        opt = *longopts;
        param = parameters;

        for ( i=0; i<nboptions; i++, opt++, param++ ) {

            /* Look for a valid option */
            while ( !(param->flags & PARAM_OPTION) ||
                    (param->name == NULL) )
            {
                param++;
            }

            opt->name    = param->name;
            opt->has_arg = ( param->has_arg > 0 ) ? 1 : 0;
            opt->flag    = NULL;
            opt->val     = param->shname;
        }
    }
#endif
}

parameter_t *
parameters_get( int shname )
{
    parameter_t *param = parameters;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            param++;
            continue;
        }

        if ( shname == param->shname ) {
            return param;
        }
        param++;
    }

    fprintf( stderr, "parameters_get could not find parameter %d(%c)\n", shname, shname );
    return NULL;
}

int
parameters_getvalue_int( const char *name )
{
    parameter_t *param = parameters;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            param++;
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            param++;
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_int should not be called with parameter %s\n", name );
            return -1;
        }

        if ( param->valtype != TestValInt ) {
            fprintf( stderr, "parameters_getvalue_int has been called with a non integer parameter (%s)\n", name );
            return -1;
        }

        return param->value.ival;
    }

    fprintf( stderr, "parameters_getvalue_int could not find parameter %s\n", name );
    return -1;
}

double
parameters_getvalue_fixdbl( const char *name )
{
    parameter_t *param = parameters;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            param++;
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            param++;
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_fixdbl should not be called with parameter %s\n", name );
            return -1;
        }

        if ( param->valtype != TestValDouble ) {
            fprintf( stderr, "parameters_getvalue_fixdbl has been called with a non float parameter (%s)\n", name );
            return -1;
        }

        return param->value.dval;
    }

    fprintf( stderr, "parameters_getvalue_fixdbl could not find parameter %s\n", name );
    return -1;
}

char *
parameters_getvalue_str( const char *name )
{
    parameter_t *param = parameters;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            param++;
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            param++;
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_str should not be called with parameter %s\n", name );
            return NULL;
        }

        if ( param->valtype != TestString ) {
            fprintf( stderr, "parameters_getvalue_str has been called with a non string parameter (%s)\n", name );
            return NULL;
        }

        return param->value.str;
    }

    fprintf( stderr, "parameters_getvalue_str could not find parameter %s\n", name );
    return NULL;
}

parameter_t *
parameters_getbyname( const char *name )
{
    parameter_t *param = parameters;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            param++;
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            param++;
            continue;
        }

        /* if ( param->has_arg < 2 ) { */
        /*     fprintf( stderr, "parameters_getbyname should not be called with parameter %s\n", name ); */
        /*     return NULL; */
        /* } */

        return param;
    }

    fprintf( stderr, "parameters_getbyname could not find parameter %s\n", name );
    return NULL;
}

void parameters_parser( int argc, char **argv )
{
    int opt;
    char optstring[STR_MAX_LENGTH];
    struct option *longopts = NULL;
    parameter_t *param;

    parameters_getopt_init( optstring, &longopts );

#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1)
#else
    while ((opt = getopt(argc, argv, optstring)) != -1)
#endif
    {
        switch(opt) {
        case 'h':
            print_usage(argv[0]);
            exit(0);

        case '?': /* error from getopt[_long] */
            exit(1);
            break;

        default:
            param = parameters_get( opt );
            if ( param == NULL ) {
                print_usage(argv[0]);
                exit(1);
            }
            parameters_addvalues( param, optarg );
        }
    }

    if ( longopts != NULL ) {
        free( longopts );
    }

#if !defined(CHAMELEON_TESTINGS_VENDOR)
    /* Force Async if splitsub is enabled */
    {
        int splitsub = parameters_getvalue_int( "splitsub" );

        if ( splitsub ) {
            param = parameters_get( 's' );
            if ( param == NULL ) {
                print_usage(argv[0]);
                exit(1);
            }
            parameters_addvalues( param, NULL );

#if defined(CHAMELEON_RUNTIME_SYNC)
            fprintf( stderr, "Spliting the submission and the execution stages is not possible when the option CHAMELEON_RUNTIME_SYNC is enabled\n" );
            exit(0);
#endif
        }
    }
#endif
}

void
parameters_destroy()
{
    parameter_t *param = parameters;
    vallist_t *current, *next;

    while( param->helper != NULL )
    {
        /* This is not an option, we skip it */
        if ( param->has_arg < 2 ) {
            param++;
            continue;
        }

        current = param->vallist;
        while ( current != NULL )
        {
            next = current->next;
            free( current );
            current = next;
        }
        param++;
    }
    return;
}
