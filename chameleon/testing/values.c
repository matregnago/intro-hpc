/**
 *
 * @file values.c
 *
 * @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon testing values toutine to read/print the parameters
 *
 * @version 1.4.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 */
#include "testings.h"
#include <math.h>

/**
 * @brief Convert the input string to an integer
 * @param[in] str
 *    The input string
 * @return The integer read.
 */
val_t pread_int( const char *str )
{
    val_t val;
    val.ival = atoi( str );
    return val;
}

/**
 * @brief Convert the input string to a float
 * @param[in] str
 *    The input string
 * @return The float read.
 */
val_t pread_float( const char *str )
{
    val_t val;
    val.sval = strtof( str, NULL );
    return val;
}

/**
 * @brief Convert the input string to a double
 * @param[in] str
 *    The input string
 * @return The double read.
 */
val_t pread_double( const char *str )
{
    val_t val;
    val.dval = strtod( str, NULL );
    return val;
}

/**
 * @brief Convert the input string to a complex single precision
 * @param[in] str
 *    The input string
 * @return The complex single precision read.
 */
val_t pread_complex32( const char *str )
{
    float re, im;
    val_t val;
    int rc;

    rc = sscanf( str, "%e,%e", &re, &im );

    if ( rc == 2 ) {
        val.cval = re + I * im;
    }
    else if (rc == 1){
        val.cval = re;
    }
    else {
        val.cval = nan("NaN");
    }

    return val;
}

/**
 * @brief Convert the input string to a complex double precision
 * @param[in] str
 *    The input string
 * @return The complex double precision read.
 */
val_t pread_complex64( const char *str )
{
    double re, im;
    val_t val;
    int rc;

    rc = sscanf( str, "%le,%le", &re, &im );

    if ( rc == 2 ) {
        val.zval = re + I * im;
    }
    else if (rc == 1){
        val.zval = re;
    }
    else {
        val.zval = nan("NaN");
    }

    return val;
}

/**
 * @brief Convert the input string to a cham_trans_t
 * @param[in] str
 *    The input string
 * @return The cham_trans_t read.
 */
val_t pread_trans( const char *str )
{
    val_t val;
    val.trans = ChamNoTrans;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamConjTrans", str ) == 0 ) ||
         ( strcasecmp( "ConjTrans", str ) == 0 ) )
    {
        val.trans = ChamConjTrans;
    }
    else if ( ( strcasecmp( "ChamTrans", str ) == 0 ) ||
              ( strcasecmp( "Trans", str ) == 0 ) )
    {
        val.trans = ChamTrans;
    }
    else if ( ( strcasecmp( "ChamNoTrans", str ) == 0 ) ||
              ( strcasecmp( "NoTrans", str ) == 0 ) )
    {
        val.trans = ChamNoTrans;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamConjTrans) || (v == (ChamConjTrans-ChamNoTrans)) ) {
            val.trans = ChamConjTrans;
        }
        else if ( (v == ChamTrans) || (v == (ChamTrans-ChamNoTrans)) ) {
            val.trans = ChamTrans;
        }
        else {
            val.trans = ChamNoTrans;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_uplo_t
 * @param[in] str
 *    The input string
 * @return The cham_uplo_t read.
 */
val_t pread_uplo( const char *str )
{
    val_t val;
    val.uplo = ChamUpperLower;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamUpper", str ) == 0 ) ||
         ( strcasecmp( "Upper",     str ) == 0 ) )
    {
        val.uplo = ChamUpper;
    }
    else if ( ( strcasecmp( "ChamLower", str ) == 0 ) ||
              ( strcasecmp( "Lower",     str ) == 0 ) )
    {
        val.uplo = ChamLower;
    }
    else if ( ( strcasecmp( "ChamUpperLower", str ) == 0 ) ||
              ( strcasecmp( "UpperLower",     str ) == 0 ) ||
              ( strcasecmp( "General",        str ) == 0 ) )
    {
        val.uplo = ChamUpperLower;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamUpper) || (v == 0) ) {
            val.uplo = ChamUpper;
        }
        else if ( (v == ChamLower) || (v == (ChamLower-ChamUpper)) ) {
            val.uplo = ChamLower;
        }
        else {
            val.uplo = ChamUpperLower;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_diag_t
 * @param[in] str
 *    The input string
 * @return The cham_diag_t read.
 */
val_t pread_diag( const char *str )
{
    val_t val;
    val.diag = ChamNonUnit;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamNonUnit", str ) == 0 ) ||
         ( strcasecmp( "NonUnit",     str ) == 0 ) )
    {
        val.diag = ChamNonUnit;
    }
    else if ( ( strcasecmp( "ChamUnit", str ) == 0 ) ||
              ( strcasecmp( "Unit",     str ) == 0 ) )
    {
        val.diag = ChamUnit;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamUnit) || (v == (ChamUnit-ChamNonUnit)) ) {
            val.diag = ChamUnit;
        }
        else {
            val.diag = ChamNonUnit;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_side_t
 * @param[in] str
 *    The input string
 * @return The cham_side_t read.
 */
val_t pread_side( const char *str )
{
    val_t val;
    val.side = ChamLeft;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamLeft", str ) == 0 ) ||
         ( strcasecmp( "Left",     str ) == 0 ) )
    {
        val.side = ChamLeft;
    }
    else if ( ( strcasecmp( "ChamRight", str ) == 0 ) ||
              ( strcasecmp( "Right",     str ) == 0 ) )
    {
        val.side = ChamRight;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamRight) || (v == (ChamRight-ChamLeft)) ) {
            val.side = ChamRight;
        }
        else {
            val.side = ChamLeft;
        }
    }
    return val;
}


/**
 * @brief Convert the input string to a cham_dir_t
 * @param[in] str
 *    The input string
 * @return The cham_dir_t read.
 */
val_t pread_dir( const char *str )
{
    val_t val;
    val.dir = ChamDirForward;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamDirForward", str ) == 0 ) ||
         ( strcasecmp( "Forward",        str ) == 0 ) )
    {
        val.dir = ChamDirForward;
    }
    else if ( ( strcasecmp( "ChamDirBackward", str ) == 0 ) ||
              ( strcasecmp( "Backward",        str ) == 0 ) )
    {
        val.dir = ChamDirBackward;
    }
    else {
        int v = atoi( str );
        if ( ( v == ChamDirForward ) || ( v == ( ChamDirForward - ChamDirBackward ) ) ) {
            val.dir = ChamDirForward;
        }
        else {
            val.dir = ChamDirBackward;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_job_t
 * @param[in] str
 *    The input string
 * @return The cham_job_t read.
 */
val_t pread_job( const char *str )
{
    val_t val;
    val.job = ChamNoVec;

    /* Keep in sync with help documentation in testing/chameleon_ztesting.c */
    if ( ( strcasecmp( "ChamNoVec", str ) == 0 ) ||
         ( strcasecmp( "NoVec",     str ) == 0 ) ||
         ( strcasecmp( "N",         str ) == 0 ) )
    {
        val.job = ChamNoVec;
    }
    else if ( ( strcasecmp( "ChamVec", str ) == 0 ) ||
              ( strcasecmp( "Vec",     str ) == 0 ) )
    {
        val.job = ChamVec;
    }
    else if ( ( strcasecmp( "ChamIvec", str ) == 0 ) ||
              ( strcasecmp( "Ivec",     str ) == 0 ) ||
              ( strcasecmp( "I",        str ) == 0 ) )
    {
        val.job = ChamIvec;
    }
    else if ( ( strcasecmp( "ChamAllVec", str ) == 0 ) ||
              ( strcasecmp( "AllVec",     str ) == 0 ) ||
              ( strcasecmp( "A",          str ) == 0 ) )
    {
        val.job = ChamAllVec;
    }
    else if ( ( strcasecmp( "ChamSVec", str ) == 0 ) ||
              ( strcasecmp( "SVec",     str ) == 0 ) ||
              ( strcasecmp( "S",        str ) == 0 ) )
    {
        val.job = ChamSVec;
    }
    else if ( ( strcasecmp( "ChamOVec", str ) == 0 ) ||
              ( strcasecmp( "OVec",     str ) == 0 ) ||
              ( strcasecmp( "O",        str ) == 0 ) )
    {
        val.job = ChamOVec;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamVec) || (v == (ChamVec-ChamNoVec)) ) {
            val.job = ChamVec;
        }
        else if ( (v == ChamIvec) || (v == (ChamIvec-ChamNoVec)) ) {
            val.job = ChamIvec;
        }
        else if ( (v == ChamAllVec) || (v == (ChamAllVec-ChamNoVec)) ) {
            val.job = ChamAllVec;
        }
        else if ( (v == ChamSVec) || (v == (ChamSVec-ChamNoVec)) ) {
            val.job = ChamSVec;
        }
        else if ( (v == ChamOVec) || (v == (ChamOVec-ChamNoVec)) ) {
            val.job = ChamOVec;
        }
        else {
            val.job = ChamNoVec;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_normtype_t
 * @param[in] str
 *    The input string
 * @return The cham_normtype_t read.
 */
val_t pread_norm( const char *str )
{
    val_t val;
    val.ntype = ChamOneNorm;

    /* Keep in sync with help documentation in testing/{chameleon,vendor}_ztesting.c */
    if ( ( strcasecmp( "ChamOneNorm", str ) == 0 ) ||
         ( strcasecmp( "OneNorm",     str ) == 0 ) ||
         ( strcasecmp( "One",         str ) == 0 ) )
    {
        val.ntype = ChamOneNorm;
    }
    else if ( ( strcasecmp( "ChamFrobeniusNorm", str ) == 0 ) ||
              ( strcasecmp( "FrobeniusNorm",     str ) == 0 ) ||
              ( strcasecmp( "Frobenius",         str ) == 0 ) ||
              ( strcasecmp( "Fro",               str ) == 0 ) )
    {
        val.ntype = ChamFrobeniusNorm;
    }
    else if ( ( strcasecmp( "ChamInfNorm", str ) == 0 ) ||
              ( strcasecmp( "InfNorm",     str ) == 0 ) ||
              ( strcasecmp( "Inf",         str ) == 0 ) )
    {
        val.ntype = ChamInfNorm;
    }
    else if ( ( strcasecmp( "ChamMaxNorm", str ) == 0 ) ||
              ( strcasecmp( "MaxNorm",     str ) == 0 ) ||
              ( strcasecmp( "Max",         str ) == 0 ) )
    {
        val.ntype = ChamMaxNorm;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamMaxNorm) || (v == (ChamMaxNorm-ChamOneNorm)) ) {
            val.ntype = ChamMaxNorm;
        }
        else if ( (v == ChamInfNorm) || (v == (ChamInfNorm-ChamOneNorm)) ) {
            val.ntype = ChamInfNorm;
        }
        else if ( (v == ChamFrobeniusNorm) || (v == (ChamFrobeniusNorm-ChamOneNorm)) ) {
            val.ntype = ChamFrobeniusNorm;
        }
        else {
            val.ntype = ChamOneNorm;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a cham_rec_t
 * @param[in] str
 *    The input string
 * @return The cham_rec_t read.
 */
val_t pread_rec( const char *str )
{
    val_t val;
    val.rec = ChamRecFull;
    if ( strcasecmp("full", str) == 0 )
    {
        val.rec = ChamRecFull;
    }
    else if ( strcasecmp("random", str) == 0 )
    {
        val.rec = ChamRecRandom;
    }
    else if ( strcasecmp("diag", str) == 0 )
    {
        val.rec = ChamRecDiag;
    }
    else if ( strcasecmp("smart", str) == 0 )
    {
        val.rec = ChamRecSmart;
    }
    else {
        int v = atoi( str );
        if ( (v == ChamRecRandom) || (v == (ChamRecRandom-ChamRecFull)) ) {
            val.rec = ChamRecRandom;
        }
        else if ( (v == ChamRecDiag) || (v == (ChamRecDiag-ChamRecFull)) ) {
            val.rec = ChamRecDiag;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
val_t pread_string( const char *str )
{
    val_t val;
    int i, len = strlen( str );
    val.str = strdup( str );
    for( i=0; i<len; i++ ) {
        if ( (val.str[i] == '\n') ||
             (val.str[i] == ' ') )
        {
            val.str[i] = '\0';
            break;
        }
    }
    return val;
}

/**
 * @brief Convert the input string to an integer
 * @param[in] str
 *    The input string
 * @return The integer read.
 */
char *sprint_int( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*d", nbchar, val.ival );
    }
    else {
        rc = sprintf( str_in, ";%d", val.ival );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a float
 * @param[in] str
 *    The input string
 * @return The float read.
 */
char *sprint_float( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*e", chameleon_max( 13, nbchar ), val.sval );
    }
    else {
        rc = sprintf( str_in, ";%e", val.sval );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a double
 * @param[in] str
 *    The input string
 * @return The double read.
 */
char *sprint_double( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*e", chameleon_max( 13, nbchar ), val.dval );
    }
    else {
        rc = sprintf( str_in, ";%e", val.dval );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a complex single precision
 * @param[in] str
 *    The input string
 * @return The complex single precision read.
 */
char *sprint_complex32( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " (%*e,%*e)", chameleon_max( 13, nbchar ), crealf(val.cval), chameleon_max( 13, nbchar ), cimagf(val.cval) );
    }
    else {
        rc = sprintf( str_in, ";%e,%e", crealf(val.cval), cimagf(val.cval) );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a complex double precision
 * @param[in] str
 *    The input string
 * @return The complex double precision read.
 */
char *sprint_complex64( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " (%*e,%*e)", chameleon_max( 13, nbchar ), creal(val.zval), chameleon_max( 13, nbchar ), cimag(val.zval) );
    }
    else {
        rc = sprintf( str_in, ";%e,%e", creal(val.zval), cimag(val.zval) );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_trans_t
 * @param[in] str
 *    The input string
 * @return The cham_trans_t read.
 */
char *sprint_trans( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.trans == ChamConjTrans) ? "ConjTrans" :
                      ((val.trans == ChamTrans) ? "Trans" : "NoTrans") );
    }
    else {
        rc = sprintf( str_in, ";%d", val.trans );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_uplo_t
 * @param[in] str
 *    The input string
 * @return The cham_uplo_t read.
 */
char *sprint_uplo( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.uplo == ChamUpper) ? "Upper" :
                      ((val.uplo == ChamLower) ? "Lower" : "General") );
    }
    else {
        rc = sprintf( str_in, ";%d", val.uplo );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_diag_t
 * @param[in] str
 *    The input string
 * @return The cham_diag_t read.
 */
char *sprint_diag( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.diag == ChamUnit) ? "Unit" : "NonUnit" );
    }
    else {
        rc = sprintf( str_in, ";%d", val.diag );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_side_t
 * @param[in] str
 *    The input string
 * @return The cham_side_t read.
 */
char *sprint_side( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.side == ChamLeft) ? "Left" : "Right" );
    }
    else {
        rc = sprintf( str_in, ";%d", val.side );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_dir_t
 * @param[in] str
 *    The input string
 * @return The cham_dir_t read.
 */
char *sprint_dir( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar,
                      (val.dir == ChamDirForward) ? "Forward" : "Backward" );
    }
    else {
        rc = sprintf( str_in, ";%d", val.dir );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_job_t
 * @param[in] str
 *    The input string
 * @return The cham_job_t read.
 */
char *sprint_job( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        char *name;
        switch( val.job ) {
        case ChamNoVec:
            name = "No Vec";
            break;
        case ChamVec:
            name = "Vec";
            break;
        case ChamIvec:
            name = "Ivec";
            break;
        case ChamAllVec:
            name = "All Vec";
            break;
        case ChamSVec:
            name = "S Vec";
            break;
        case ChamOVec:
            name = "O Vec";
            break;
        default:
            name = "ERR";
        }
        rc = sprintf( str_in, " %-*s", nbchar, name );
    }
    else {
        rc = sprintf( str_in, ";%d", val.job );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_normtype_t
 * @param[in] str
 *    The input string
 * @return The cham_normtype_t read.
 */
char *sprint_norm( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        char *name;
        switch( val.ntype ) {
        case ChamMaxNorm:
            name = "Max";
            break;
        case ChamOneNorm:
            name = "One";
            break;
        case ChamInfNorm:
            name = "Inf";
            break;
        case ChamFrobeniusNorm:
            name = "Frb";
            break;
        default:
            name = "ERR";
        }
        rc = sprintf( str_in, " %-*s", nbchar, name );
    }
    else {
        rc = sprintf( str_in, ";%d", val.ntype );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a cham_rec_t
 * @param[in] str
 *    The input string
 * @return The cham_rec_t read.
 */
char *sprint_rec( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        char *name;
        switch( val.rec ) {
        case ChamRecFull:
            name = "Full";
            break;
        case ChamRecRandom:
            name = "Rand";
            break;
        case ChamRecDiag:
            name = "Diag";
            break;
        default:
            name = "ERR";
        }
        rc = sprintf( str_in, " %-*s", nbchar, name );
    }
    else {
        rc = sprintf( str_in, ";%d", val.ntype );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
char *sprint_check( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %*s", nbchar, ( val.ival == 0 ) ? "SUCCESS" : "FAILED" );
    }
    else {
        rc = sprintf( str_in, ";%s", ( val.ival == 0 ) ? "SUCCESS" : "FAILED" );
    }
    return str_in+rc;
}

/**
 * @brief Convert the input string to a string
 * @param[in] str
 *    The input string
 * @return The string read.
 */
char *sprint_string( val_t val, int human, int nbchar, char *str_in )
{
    int rc;
    if ( human ) {
        rc = sprintf( str_in, " %-*s", nbchar, val.str );
    }
    else {
        rc = sprintf( str_in, ";%s", val.str );
    }
    return str_in+rc;
}

#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f

/**
 * @brief Generate a random number
 */
static inline unsigned long long int
testing_Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, ++i) {
    if (n & 1)
      ran = a_k * ran + c_k;
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

/**
 * @brief Generate a random float
 */
float
testing_salea()
{
    float val;
    unsigned long long int ran;

    ran = testing_Rnd64_jump( 2, random() );

    /* Real part */
    val = 0.5f - ran * RndF_Mul;

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &val, 1, MPI_FLOAT, 0, MPI_COMM_WORLD );
#endif

    return val;
}

/**
 * @brief Generate a random integer
 */
long
testing_ialea()
{
    long r = random();
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &r, 1, MPI_LONG, 0, MPI_COMM_WORLD );
#endif
    return r;
}

/**
 * @brief Generate a random double
 */
double
testing_dalea()
{
    double val;
    unsigned long long int ran;

    ran = testing_Rnd64_jump( 2, random() );

    /* Real part */
    val = 0.5f - ran * RndF_Mul;

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

    return val;
}

/**
 * @brief Generate a random complex single precision
 */
CHAMELEON_Complex32_t
testing_calea()
{
    CHAMELEON_Complex32_t val;
    unsigned long long int ran;

    ran = testing_Rnd64_jump( 2, random() );

    /* Real part */
    val = 0.5f - ran * RndF_Mul;
    ran  = Rnd64_A * ran + Rnd64_C;

    /* Imaginary part */
    val += I*(0.5f - ran * RndF_Mul);

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &val, 1, MPI_COMPLEX, 0, MPI_COMM_WORLD );
#endif

    return val;
}

/**
 * @brief Generate a random complex double precision
 */
CHAMELEON_Complex64_t
testing_zalea()
{
    CHAMELEON_Complex64_t val;
    unsigned long long int ran;

    ran = testing_Rnd64_jump( 2, random() );

    /* Real part */
    val = 0.5f - ran * RndF_Mul;
    ran  = Rnd64_A * ran + Rnd64_C;

    /* Imaginary part */
    val += I*(0.5f - ran * RndF_Mul);

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast( &val, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );
#endif

    return val;
}
