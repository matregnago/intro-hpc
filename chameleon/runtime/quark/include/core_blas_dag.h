/**
 *
 * @file quark/core_blas_dag.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark DAG generation header
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 *
 */
#ifndef _core_blas_dag_h_
#define _core_blas_dag_h_

#if defined(QUARK_DOT_DAG_ENABLE) /* || 1 */
#define DAG_SET_PROPERTIES( _name, _color )                             \
    QUARK_Task_Flag_Set( (Quark_Task_Flags*)opt, TASK_LABEL, (intptr_t)(_name)  ); \
    QUARK_Task_Flag_Set( (Quark_Task_Flags*)opt, TASK_COLOR, (intptr_t)(_color) );
#else
#define DAG_SET_PROPERTIES( _name, _color )
#endif

#define DAG_CORE_ASUM       DAG_SET_PROPERTIES( "ASUM"      , "white"   )
#define DAG_CORE_AXPY       DAG_SET_PROPERTIES( "AXPY"      , "white"   )
#define DAG_CORE_BUILD      DAG_SET_PROPERTIES( "BUILD"     , "white"   )
#define DAG_CORE_GEADD      DAG_SET_PROPERTIES( "GEADD"     , "white"   )
#define DAG_CORE_LASCAL     DAG_SET_PROPERTIES( "LASCAL"    , "white"   )
#define DAG_CORE_GELQT      DAG_SET_PROPERTIES( "GELQT"     , "green"   )
#define DAG_CORE_GEMM       DAG_SET_PROPERTIES( "GEMM"      , "yellow"  )
#define DAG_CORE_GEQRT      DAG_SET_PROPERTIES( "GEQRT"     , "green"   )
#define DAG_CORE_GESSM      DAG_SET_PROPERTIES( "GESSM"     , "cyan"    )
#define DAG_CORE_GETRF      DAG_SET_PROPERTIES( "GETRF"     , "green"   )
#define DAG_CORE_GETRIP     DAG_SET_PROPERTIES( "GETRIP"    , "white"   )
#define DAG_CORE_LATRO      DAG_SET_PROPERTIES( "LATRO"     , "white"   )
#define DAG_CORE_HEMM       DAG_SET_PROPERTIES( "HEMM"      , "white"   )
#define DAG_CORE_HER2K      DAG_SET_PROPERTIES( "HER2K"     , "white"   )
#define DAG_CORE_HERFB      DAG_SET_PROPERTIES( "HERFB"     , "pink"    )
#define DAG_CORE_HERK       DAG_SET_PROPERTIES( "HERK"      , "yellow"  )
#define DAG_CORE_LACPY      DAG_SET_PROPERTIES( "LACPY"     , "white"   )
#define DAG_CORE_LAG2C      DAG_SET_PROPERTIES( "LAG2C"     , "white"   )
#define DAG_CORE_LAG2Z      DAG_SET_PROPERTIES( "LAG2Z"     , "white"   )
#define DAG_CORE_LANGE      DAG_SET_PROPERTIES( "LANGE"     , "white"   )
#define DAG_CORE_LANGE_MAX  DAG_SET_PROPERTIES( "LANGE_MAX" , "white"   )
#define DAG_CORE_LANHE      DAG_SET_PROPERTIES( "LANHE"     , "white"   )
#define DAG_CORE_LANSY      DAG_SET_PROPERTIES( "LANSY"     , "white"   )
#define DAG_CORE_LANTR      DAG_SET_PROPERTIES( "LANTR"     , "white"   )
#define DAG_CORE_LASET      DAG_SET_PROPERTIES( "LASET"     , "orange"  )
#define DAG_CORE_LASWP      DAG_SET_PROPERTIES( "LASWP"     , "orange"  )
#define DAG_CORE_LAUUM      DAG_SET_PROPERTIES( "LAUUM"     , "white"   )
#define DAG_CORE_PLGHE      DAG_SET_PROPERTIES( "PLGHE"     , "white"   )
#define DAG_CORE_PLGSY      DAG_SET_PROPERTIES( "PLGSY"     , "white"   )
#define DAG_CORE_PLRNT      DAG_SET_PROPERTIES( "PLRNT"     , "white"   )
#define DAG_CORE_LASSQ      DAG_SET_PROPERTIES( "LASSQ"     , "white"   )
#define DAG_CORE_POTRF      DAG_SET_PROPERTIES( "POTRF"     , "green"   )
#define DAG_CORE_SHIFT      DAG_SET_PROPERTIES( "SHIFT"     , "white"   )
#define DAG_CORE_SHIFTW     DAG_SET_PROPERTIES( "SHIFTW"    , "white"   )
#define DAG_CORE_SSSSM      DAG_SET_PROPERTIES( "SSSSM"     , "yellow"  )
#define DAG_CORE_SWPAB      DAG_SET_PROPERTIES( "SWPAB"     , "white"   )
#define DAG_CORE_SYMM       DAG_SET_PROPERTIES( "SYMM"      , "white"   )
#define DAG_CORE_SYR2K      DAG_SET_PROPERTIES( "SYR2K"     , "white"   )
#define DAG_CORE_SYRK       DAG_SET_PROPERTIES( "SYRK"      , "red"     )
#define DAG_CORE_TRMM       DAG_SET_PROPERTIES( "TRMM"      , "cyan"    )
#define DAG_CORE_TRSM       DAG_SET_PROPERTIES( "TRSM"      , "cyan"    )
#define DAG_CORE_TRSSQ      DAG_SET_PROPERTIES( "TRSSQ"     , "white"   )
#define DAG_CORE_TRTRI      DAG_SET_PROPERTIES( "TRTRI"     , "white"   )
#define DAG_CORE_TPLQT      DAG_SET_PROPERTIES( "TPLQT"     , "red"     )
#define DAG_CORE_TPMLQT     DAG_SET_PROPERTIES( "TPMLQT"    , "yellow"  )
#define DAG_CORE_TPMQRT     DAG_SET_PROPERTIES( "TPMQRT"    , "yellow"  )
#define DAG_CORE_TPQRT      DAG_SET_PROPERTIES( "TPQRT"     , "red"     )
#define DAG_CORE_TSTRF      DAG_SET_PROPERTIES( "TSTRF"     , "red"     )
#define DAG_CORE_UNMLQ      DAG_SET_PROPERTIES( "UNMLQ"     , "cyan"    )
#define DAG_CORE_UNMQR      DAG_SET_PROPERTIES( "UNMQR"     , "cyan"    )
#define DAG_CORE_ORMLQ      DAG_SET_PROPERTIES( "ORMLQ"     , "cyan"    )
#define DAG_CORE_ORMQR      DAG_SET_PROPERTIES( "ORMQR"     , "cyan"    )

#define DAG_CORE_TSLQT      DAG_CORE_TPLQT
#define DAG_CORE_TSMLQ      DAG_CORE_TPMLQT
#define DAG_CORE_TSMQR      DAG_CORE_TPMQRT
#define DAG_CORE_TSQRT      DAG_CORE_TPQRT
#define DAG_CORE_TTLQT      DAG_CORE_TPLQT
#define DAG_CORE_TTMLQ      DAG_CORE_TPMLQT
#define DAG_CORE_TTMQR      DAG_CORE_TPMQRT
#define DAG_CORE_TTQRT      DAG_CORE_TPQRT

#define DAG_CORE_GESSQ      DAG_SET_PROPERTIES( "GESSQ"    , "white"    )
#define DAG_CORE_HESSQ      DAG_SET_PROPERTIES( "HESSQ"    , "white"    )
#define DAG_CORE_SYSSQ      DAG_SET_PROPERTIES( "SYSSQ"    , "white"    )
#define DAG_CORE_TRSSQ      DAG_SET_PROPERTIES( "TRSSQ"    , "white"    )
#define DAG_CORE_PLSSQ      DAG_SET_PROPERTIES( "PLSSQ"    , "white"    )
#define DAG_CORE_PLSSQ2     DAG_SET_PROPERTIES( "PLSSQ2"   , "white"    )

#define DAG_CORE_GESUM      DAG_SET_PROPERTIES( "GESUM"    , "white"    )
#define DAG_CORE_CESCA      DAG_SET_PROPERTIES( "CESCA"    , "orange"   )
#define DAG_CORE_GRAM       DAG_SET_PROPERTIES( "GRAM"     , "orange"   )

#endif /* _core_blas_dag_h_ */
