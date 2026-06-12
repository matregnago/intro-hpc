/**
 *
 * @file chameleon_constants.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global constants
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @author Terry Cojean
 * @author Matthieu Kuhn
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_constants_h_
#define _chameleon_constants_h_

/**
 *
 * @brief Chameleon constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 */
/**
 * @brief Matrix floating point arithmetic
 */
typedef enum chameleon_arithmetic_e {
    Cham8      = 0,
    ChamHalf   = 1,
    ChamSingle = 2,
    ChamDouble = 3,
} cham_arithmetic_t;

typedef enum chameleon_ftype_e {
    ChamInt     = 0,
    ChamReal    = 1,
    ChamComplex = 2,
} cham_ftype_t;

/**
 * Defines the arithmetic type mask using the 1st and 2nd bits
 */
#define CHAM_ARITHMETIC_MASK ( 0x3 )
/**
 * Defines the floating type mask using the 3rd and 4th bits
 */
#define CHAM_FTYPE_MASK ( 0x3 << 2 )
/**
 * Defines the mixed precision flag mask on the 5th bit
 */
#define CHAM_MIXED_MASK ( 0x1 << 4 )

#define cham_get_arith( _ftype_ )   ( (_ftype_) & CHAM_ARITHMETIC_MASK )
#define cham_get_ftype( _ftype_ )   (( (_ftype_) & CHAM_FTYPE_MASK ) >> 2 )
#define cham_get_flttype( _ftype_ ) ( (_ftype_) & (CHAM_FTYPE_MASK | CHAM_ARITHMETIC_MASK) )
#define cham_is_mixed( _ftype_ )    ( (_ftype_) & CHAM_MIXED_MASK )

#define cham_clean_mixed( _ftype_ )    ( (_ftype_) & ~CHAM_MIXED_MASK )

#define CHAMELEON_FLTTYPE( _ftype_, _arithmetic_ ) ( ((_ftype_) << 2) | (_arithmetic_) )

typedef enum chameleon_flttype_e {
    ChamByte               = CHAMELEON_FLTTYPE( ChamInt,     Cham8      ),
    ChamInteger16          = CHAMELEON_FLTTYPE( ChamInt,     ChamHalf   ),
    ChamInteger            = CHAMELEON_FLTTYPE( ChamInt,     ChamSingle ),
    ChamInteger32          = CHAMELEON_FLTTYPE( ChamInt,     ChamSingle ),
    ChamInteger64          = CHAMELEON_FLTTYPE( ChamInt,     ChamDouble ),
    ChamRealHalf           = CHAMELEON_FLTTYPE( ChamReal,    ChamHalf   ),
    ChamRealFloat          = CHAMELEON_FLTTYPE( ChamReal,    ChamSingle ),
    ChamRealDouble         = CHAMELEON_FLTTYPE( ChamReal,    ChamDouble ),
    ChamComplexHalf        = CHAMELEON_FLTTYPE( ChamComplex, ChamHalf   ),
    ChamComplexFloat       = CHAMELEON_FLTTYPE( ChamComplex, ChamSingle ),
    ChamComplexDouble      = CHAMELEON_FLTTYPE( ChamComplex, ChamDouble ),
    ChamRealDoubleMixed    = ChamRealDouble    | CHAM_MIXED_MASK,
    ChamComplexDoubleMixed = ChamComplexDouble | CHAM_MIXED_MASK,
} cham_flttype_t;

#define ChamComplexSingle ChamComplexFloat
#define ChamRealSingle    ChamRealFloat

#define ChamConvert( in, out ) ( cham_clean_mixed(in) | (cham_clean_mixed(out) << 5) )

#define ChamConvertComplexDoubleToDouble ChamConvert( ChamComplexDouble, ChamComplexDouble )
#define ChamConvertComplexDoubleToSingle ChamConvert( ChamComplexDouble, ChamComplexSingle )
#define ChamConvertComplexDoubleToHalf   ChamConvert( ChamComplexDouble, ChamComplexHalf   )
#define ChamConvertComplexSingleToDouble ChamConvert( ChamComplexSingle, ChamComplexDouble )
#define ChamConvertComplexSingleToSingle ChamConvert( ChamComplexSingle, ChamComplexSingle )
#define ChamConvertComplexSingleToHalf   ChamConvert( ChamComplexSingle, ChamComplexHalf   )
#define ChamConvertComplexHalfToDouble   ChamConvert( ChamComplexHalf,   ChamComplexDouble )
#define ChamConvertComplexHalfToSingle   ChamConvert( ChamComplexHalf,   ChamComplexSingle )
#define ChamConvertComplexHalfToHalf     ChamConvert( ChamComplexHalf,   ChamComplexHalf   )

#define ChamConvertRealDoubleToDouble ChamConvert( ChamRealDouble, ChamRealDouble )
#define ChamConvertRealDoubleToSingle ChamConvert( ChamRealDouble, ChamRealSingle )
#define ChamConvertRealDoubleToHalf   ChamConvert( ChamRealDouble, ChamRealHalf   )
#define ChamConvertRealSingleToDouble ChamConvert( ChamRealSingle, ChamRealDouble )
#define ChamConvertRealSingleToSingle ChamConvert( ChamRealSingle, ChamRealSingle )
#define ChamConvertRealSingleToHalf   ChamConvert( ChamRealSingle, ChamRealHalf   )
#define ChamConvertRealHalfToDouble   ChamConvert( ChamRealHalf,   ChamRealDouble )
#define ChamConvertRealHalfToSingle   ChamConvert( ChamRealHalf,   ChamRealSingle )
#define ChamConvertRealHalfToHalf     ChamConvert( ChamRealHalf,   ChamRealHalf   )

/**
 * @brief Matrix dimensions naming
 */
typedef enum chameleon_dim_e {
    DIM_m = 0,
    DIM_n = 1,
} cham_dim_t;

/**
 * @brief Matrix tile storage
 */
typedef enum chameleon_storage_e {
    ChamCM   = 101,
    ChamRM   = 102,
    ChamCCRB = 103,
    ChamCRRB = 104,
    ChamRCRB = 105,
    ChamRRRB = 106,
} cham_storage_t;

/**
 * @brief Transpostion
 */
typedef enum chameleon_trans_e {
    ChamNoTrans     = 111, /**< Use A         */
    ChamTrans       = 112, /**< Use A^t       */
    ChamConjTrans   = 113, /**< Use conj(A^t) */
    ChamConjNoTrans = 114, /**< Use conj(A)   */
} cham_trans_t;

static inline int
isValidTrans( cham_trans_t trans )
{
    if ( (trans >= ChamNoTrans) && (trans <= ChamConjNoTrans) ) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief Upper/Lower part
 */
typedef enum chameleon_uplo_e {
    ChamUpper      = 121, /**< Use lower triangle of A */
    ChamLower      = 122, /**< Use upper triangle of A */
    ChamUpperLower = 123  /**< Use the full A          */
} cham_uplo_t;

/**
 * @brief Diagonal
 */
typedef enum chameleon_diag_e {
    ChamNonUnit = 131, /**< Diagonal is non unitary */
    ChamUnit    = 132  /**< Diagonal is unitary     */
} cham_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum chameleon_side_e {
    ChamLeft  = 141, /**< Apply operator on the left  */
    ChamRight = 142  /**< Apply operator on the right */
} cham_side_t;

/**
 * @brief Norms
 */
typedef enum chameleon_normtype_e {
    ChamOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    /* ChamRealOneNorm   = 172, */
    ChamTwoNorm       = 173, /**< Two norm: max( sigma_i )                     */
    ChamFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    ChamInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    /* ChamRealInfNorm   = 176, */
    ChamMaxNorm       = 177, /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
    /* ChamRealMaxNorm   = 178 */
} cham_normtype_t;

/**
 * @brief Recursive algorithm used to partition the matrix
 */
typedef enum chameleon_rec_e {
    ChamRecFull   = 181, /**< Fully partitioned matrix */
    ChamRecRandom = 182, /**< Randomly partitioned matrix */
    ChamRecDiag   = 183, /**< Diagonaly partitioned matrix */
    ChamRecSmart  = 184, /**< Partitioning optimized for GEMM */
} cham_rec_t;

/**
 * @brief Random distribution for matrix generator
 */
typedef enum chameleon_dist_e {
    ChamDistUniform   = 201,
    ChamDistSymmetric = 202,
    ChamDistNormal    = 203,
} cham_dist_t;

/**
 * @brief Main matrix type
 */
typedef enum chameleon_mtxtype_e {
    ChamGeneral    = 231,
    ChamSymmetric  = 232,
    ChamHermitian  = 233,
    ChamTriangular = 234,
} cham_mtxtype_t;

/**
 * @brief Eigen and singular values generator format
 */
typedef enum chameleon_sym_e {
    ChamHermGeev   = 241,
    ChamHermPoev   = 242,
    ChamNonsymPosv = 243,
    ChamSymPosv    = 244
} cham_sym_t;

#define ChamNoPacking       291
#define ChamPackSubdiag     292
#define ChamPackSupdiag     293
#define ChamPackColumn      294
#define ChamPackRow         295
#define ChamPackLowerBand   296
#define ChamPackUpeprBand   297
#define ChamPackAll         298

/**
 * @brief Singular/Eigen vector job description
 */
typedef enum chameleon_job_e {
    ChamNoVec  = 301,
    ChamVec    = 302,
    ChamIvec   = 303,
    ChamAllVec = 304,
    ChamSVec   = 305,
    ChamOVec   = 306,
} cham_job_t;

/**
 * @brief Algorithm Direction
 */
typedef enum chameleon_dir_e {
    ChamDirForward  = 391, /**< Forward direction   */
    ChamDirBackward = 392, /**< Backward direction  */
} cham_dir_t;

/**
 * @brief Direction of the main vectors as for the householder reflectors in QR/LQ factorizations.
 */
typedef enum chameleon_store_e {
    ChamColumnwise = 401, /**< Column wise storage  */
    ChamRowwise    = 402, /**< Row wise storage     */
    ChamEltwise    = 403, /**< Element by element storage */
} cham_store_t;

/**
 * @brief Data access types.
 */
typedef enum chameleon_access_e {
    ChamR  = (1 << 0),        /**< Read only  */
    ChamW  = (1 << 1),        /**< Write only */
    ChamRW = (ChamR | ChamW), /**< Read-Write */
} cham_access_t;

/**
 * @brief Chameleon Bcast algorithms
 */
typedef enum chameleon_bcast_e {
    ChamBcastFull = 0, /**< Broadcast basic using a flat tree        */
    ChamBcastRing = 1, /**< Broadcast using a ring of communications */
    ChamBcastNumber    /**< Number of Bcast algorithms               */
} cham_bcast_t;

/**
 * @brief Chameleon GEMM-like algorithms
 */
typedef enum chameleon_gemm_e {
    ChamGemmAlgAuto    = -1,
    ChamGemmAlgGeneric =  0,
    ChamGemmAlgSummaA  =  1,
    ChamGemmAlgSummaB  =  2,
    ChamGemmAlgSummaC  =  3,
} cham_gemm_t;

/**
 * @brief Chameleon GETRF algorithm variants
 */
typedef enum chameleon_getrf_e {
    ChamGetrfNoPiv          = 0,
    ChamGetrfNoPivPerColumn = 1,
    ChamGetrfPPiv           = 2,
    ChamGetrfPPivPerColumn  = 3,
} cham_getrf_t;

/**
 * @brief Chameleon GETRF all reduce algorithm variants
 */
typedef enum chameleon_getrf_allreduce_e {
    ChamStarPUTasks = 0,
} cham_getrf_allreduce_t;

#define ChameleonTrd            1001
#define ChameleonBrd            1002

#define ChameleonW               501
#define ChameleonA2              502

/**
 *  CHAMELEON constants - boolean
 */
#define CHAMELEON_FALSE  0
#define CHAMELEON_TRUE   1

#define CHAMELEON_CPU    ((1ULL)<<1)
#define CHAMELEON_CUDA   ((1ULL)<<3)

/**
 *  State machine switches
 */
#define CHAMELEON_WARNINGS            1
#define CHAMELEON_ERRORS              2
#define CHAMELEON_AUTOTUNING          3
#define CHAMELEON_DAG                 4
#define CHAMELEON_GENERATE_TRACE      5
#define CHAMELEON_PROFILING_MODE      CHAMELEON_GENERATE_TRACE  /* _deprecated_ */
#define CHAMELEON_GENERATE_STATS      6
#define CHAMELEON_KERNELPROFILE_MODE  CHAMELEON_GENERATE_STATS  /* _deprecated_ */
#define CHAMELEON_PARALLEL_KERNEL     7
#define CHAMELEON_BOUND               8
#define CHAMELEON_PROGRESS            9
#define CHAMELEON_GEMM3M             10
#define CHAMELEON_GENERIC            11

/**
 *  CHAMELEON constants - configuration parameters
 */
#define CHAMELEON_CONCURRENCY       1
#define CHAMELEON_TILE_SIZE         2
#define CHAMELEON_INNER_BLOCK_SIZE  3
#define CHAMELEON_HOUSEHOLDER_MODE  5
#define CHAMELEON_HOUSEHOLDER_SIZE  6
#define CHAMELEON_TRANSLATION_MODE  7
#define CHAMELEON_LOOKAHEAD         8
#define CHAMELEON_RUNTIME           9

/**
 * @brief QR/LQ factorization trees
 */
typedef enum chameleon_householder_e {
    ChamFlatHouseholder = 1,
    ChamTreeHouseholder = 2,
} cham_householder_t;

/**
 * @brief Translation types
 */
typedef enum chameleon_translation_e {
    ChamInPlace    = 1,
    ChamOutOfPlace = 2,
} cham_translation_t;

/**
 *  @brief Descriptor access types used in data conversion function as well as data flush
 */
typedef enum chameleon_desc_acccess_e {
    ChamDescInput  = 1 << 0,
    ChamDescOutput = 1 << 1,
    ChamDescInout  = (ChamDescInput | ChamDescOutput),
} cham_desc_access_t;

/**
 * @brief Constant to describe how to initialize the mat pointer in descriptors
 */
#define CHAMELEON_MAT_ALLOC_GLOBAL NULL
#define CHAMELEON_MAT_ALLOC_TILE   ((void*)-1)
#define CHAMELEON_MAT_OOC          ((void*)-2)

/**
 * @brief Integer constant matching the previous ones for case values
 */
#define CHAMELEON_MAT_CASE_ALLOC_GLOBAL  0
#define CHAMELEON_MAT_CASE_ALLOC_TILE   -1
#define CHAMELEON_MAT_CASE_OOC          -2

/**
 *  CHAMELEON constants - success & error codes
 */
#define CHAMELEON_SUCCESS                 0
#define CHAMELEON_ERR_NOT_INITIALIZED  -101
#define CHAMELEON_ERR_REINITIALIZED    -102
#define CHAMELEON_ERR_NOT_SUPPORTED    -103
#define CHAMELEON_ERR_ILLEGAL_VALUE    -104
#define CHAMELEON_ERR_NOT_FOUND        -105
#define CHAMELEON_ERR_OUT_OF_RESOURCES -106
#define CHAMELEON_ERR_INTERNAL_LIMIT   -107
#define CHAMELEON_ERR_UNALLOCATED      -108
#define CHAMELEON_ERR_FILESYSTEM       -109
#define CHAMELEON_ERR_UNEXPECTED       -110
#define CHAMELEON_ERR_SEQUENCE_FLUSHED -111

#endif /* _chameleon_constants_h_ */
