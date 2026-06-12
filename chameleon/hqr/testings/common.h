/**
 *
 * @file common.h
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 0.1.0
 * @author Mathieu Faverge
 * @date 2017-04-27
 *
 */
#ifndef _hqr_testings_common_h_
#define _hqr_testings_common_h_

/**
 * @brief List of integer options
 */
typedef enum iparam_e {
  IPARAM_NNODES,       /**< Number of nodes                             */
  IPARAM_NCORES,       /**< Number of cores                             */
  IPARAM_P,            /**< Rows in the process grid                    */
  IPARAM_Q,            /**< Columns in the process grid                 */
  IPARAM_MT,           /**< Number of tile rows of the matrix           */
  IPARAM_NT,           /**< Number of tile columns of the matrix        */
  IPARAM_CHECK,        /**< Checking activated or not                   */
  IPARAM_VERBOSE,      /**< How much noise do we want?                  */
  IPARAM_LOWLVL_TREE,  /**< Tree used for reduction inside nodes        */
  IPARAM_HIGHLVL_TREE, /**< Tree used for reduction between nodes       */
  IPARAM_QR_TS_SZE,    /**< Size of TS domain                           */
  IPARAM_QR_HLVL_SZE,  /**< Size of the high level tree                 */
  IPARAM_QR_DOMINO,    /**< Enable/disable the domino tree              */
  IPARAM_QR_TSRR,      /**< Enable/disable the round-robin on TS domain */
  IPARAM_SIZEOF        /**< Size of the parameter array                 */
} iparam_e;

/**
 * @brief Macro to copy/paste all testings parameters
 */
#define PASTE_CODE_IPARAM_LOCALS(iparam)                                \
    int nodes = iparam[IPARAM_NNODES];                                  \
    int cores = iparam[IPARAM_NCORES];                                  \
    int P     = iparam[IPARAM_P];                                       \
    int Q     = iparam[IPARAM_Q];                                       \
    int MT    = iparam[IPARAM_MT];                                      \
    int NT    = iparam[IPARAM_NT];                                      \
    int check = iparam[IPARAM_CHECK];                                   \
    int loud  = iparam[IPARAM_VERBOSE];                                 \
    int llvl  = iparam[IPARAM_LOWLVL_TREE];                             \
    int hlvl  = iparam[IPARAM_HIGHLVL_TREE];                            \
    int qr_a  = iparam[IPARAM_QR_TS_SZE];                               \
    int qr_p  = iparam[IPARAM_QR_HLVL_SZE];                             \
    int domino= iparam[IPARAM_QR_DOMINO];                               \
    int tsrr  = iparam[IPARAM_QR_TSRR];                                 \
    (void)nodes;                                                        \
    (void)cores;                                                        \
    (void)P;                                                            \
    (void)Q;                                                            \
    (void)MT;                                                           \
    (void)NT;                                                           \
    (void)llvl;                                                         \
    (void)hlvl;                                                         \
    (void)qr_a;                                                         \
    (void)qr_p;                                                         \
    (void)domino;                                                       \
    (void)tsrr;                                                         \
    (void)loud;                                                         \
    (void)check;

/**
 * @brief Function to parse ptions argument of the testings
 */
void parse_arguments(int *_argc, char ***_argv, int *iparam);

#endif /* _hqr_testings_common_h_ */
