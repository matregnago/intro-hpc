/**
 *
 * @file libhqr_internal.h
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#ifndef _libhqr_internal_h_
#define _libhqr_internal_h_

#include "libhqr.h"
#include <stdio.h>
#include <assert.h>

#define PRINT_PIVGEN 0
#ifdef PRINT_PIVGEN
#define myassert( test ) do { if ( ! (test) ) { return -1; } } while(0)
#else
#define myassert(test)   do { assert((test)); return -1;     } while(0)
#endif

#if defined(LIBHQR_HAVE_FALLTHROUGH)
#define libhqr_attr_fallthrough __attribute__((fallthrough))
#else
#define libhqr_attr_fallthrough do {} while(0)
#endif

/**
 * @brief Enum of the different type of tree functions
 */
typedef enum qrtree_type_ {
    LIBHQR_QRTREE_UNSET = 0, /**< Undefined QR tree */
    LIBHQR_QRTREE_HQR,       /**< Hierarchical tree for QR factorization */
    LIBHQR_QRTREE_SVD,       /**< Hierarchical tree for SVD computation  */
    LIBHQR_QRTREE_SYS,       /**< Systolic tree for QR factorization     */
    LIBHQR_QRTREE_MTX,       /**< Any kind of tree stored in a matrix    */
} qrtree_type_e;

struct hqr_args_s;
typedef struct hqr_args_s hqr_args_t;

struct hqr_subpiv_s;
typedef struct hqr_subpiv_s hqr_subpiv_t;

/**
 * @brief Argument structure to store the high and low level trees used in the hierarchical reduction tree
 */
struct hqr_args_s {
    hqr_subpiv_t *llvl; /**< Pointer to the low level tree data structure (shared memory)       */
    hqr_subpiv_t *hlvl; /**< Pointer to the high level tree data structure (distributed memory) */
};

struct hqr_subpiv_s {
    /**
     * currpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] m   line you want to eliminate
     *
     *  @return the annihilator p used with m at step k
     */
    int (*currpiv)(const hqr_subpiv_t *arg, int k, int m);
    /**
     * nextpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] p   line currently used as an annihilator
     *    @param[in] m   line actually annihilated.
     *          m = MT to find the first time p is used as an annihilator during step k
     *
     *  @return the next line m' using the line p as annihilator during step k
     *          mt if p will never be used again as an annihilator.
     */
    int (*nextpiv)(const hqr_subpiv_t *arg, int k, int p, int m);
    /**
     * prevpiv
     *    @param[in] arg pointer to the qr_piv structure
     *    @param[in] k   step in the factorization
     *    @param[in] p   line currently used as an annihilator
     *    @param[in] m   line actually annihilated.
     *          m = p to find the last time p has been used as an annihilator during step k
     *
     *  @return the previous line m' using the line p as annihilator during step k
     *          mt if p has never been used before that as an annihilator.
     */
    int (*prevpiv)(const hqr_subpiv_t *arg, int k, int p, int m);
    int *ipiv;  /**< Pivot array for tree that uses pre-computed information            */
    int minMN;  /**< Minimum tile numbers min(A->mt, A->nt)                             */
    int ldd;    /**< Leading dimension of the local subtree                             */
    int a;      /**< a parameter of the subtree to define the set of tiles killes by TS */
    int p;      /**< p parameter defining the number od distributed processes           */
    int domino; /**< domino switch to enable the domino tree in distributed             */
};

/**
 * @brief Minimal structure to store the information related to each tile.
 */
typedef struct libhqr_tile_info_s {
    int type;          /**< The type of the reduction applied to the tile (@sa libhqr_type_e) */
    int index;         /**< The index of the killer in the column if not killed by TS         */
    int currpiv;       /**< The row index of the pivot for this tile                          */
    int nextpiv;       /**< The next tile for which the currpiv is a pivot                    */
    int prevpiv;       /**< The previous tile for which currpiv was a pivot                   */
    int first_nextpiv; /**< If type != 0, the first tile for which this tile is a pivot       */
    int first_prevpiv; /**< If type != 0, the last tile for which this tile is a pivot        */
} libhqr_tile_info_t;

/**
 * @brief Arguments structure for the rdmtx driver
 */
typedef struct libhqr_tile_args_s {
    libhqr_tile_info_t *tileinfo; /**< Matrix of tile information to regenerate the reduction tree        */
    int                *pivots;   /**< Main pivot at each step of the factorization (k at step k if NULL) */
    int                *nbgeqrf;  /**< Nbr of geqrf per column                                            */
    int                *killers;  /**< Matrix array to store the indices of the geqrf at each
                                       step when it is not possible to compute them                       */
} libhqr_tile_args_t;

/**
 * @brief Return the minimum of two integers
 */
static inline int libhqr_imin(int a, int b){
    return (a > b) ? b : a;
}

/**
 * @brief Return the maximum of two integers
 */
static inline int libhqr_imax(int a, int b){
    return (a > b) ? a : b;
}

/**
 * @brief Return ceil division of two integers
 */
static inline int libhqr_iceil(int a, int b){
    return (a + b - 1) / b;
}

/**
 * @brief Number of extra tile of type 1 between the tile of type 3 and the first of nb11
 */
#define nbextra1_formula (( (k % pa) > (pa - p) ) ? (-k)%pa + pa : 0)

/**
 * @name Drawing functions
 * @{
 */
void draw_rectangle(int k, int p, int m, int step_m, FILE *file);
void draw_lines(const libhqr_tree_t *qrtree, int k, int *tiles, int *step, FILE *file);
void draw_horizontal_line(int k, int p, int m, int step_p, int step_m, FILE *file);
void draw_vertical_line(  int k, int p, int m,             int step_m, FILE *file);

/**
 * @}
 *
 * @name Low level trees
 * @{
 */
void hqr_low_flat_init     ( hqr_subpiv_t *arg );
void hqr_low_binary_init   ( hqr_subpiv_t *arg );
void hqr_low_fibonacci_init( hqr_subpiv_t *arg, int minMN );
void hqr_low_greedy_init   ( hqr_subpiv_t *arg, int minMN );
void hqr_low_greedy1p_init ( hqr_subpiv_t *arg, int minMN );
void svd_low_adaptiv_init  ( hqr_subpiv_t *arg, int gmt, int gnt, int nbcores, int ratio );

/**
 * @}
 *
 * @name High level trees
 * @{
 */
void hqr_high_flat_init     ( hqr_subpiv_t *arg );
void hqr_high_binary_init   ( hqr_subpiv_t *arg );
void hqr_high_fibonacci_init( hqr_subpiv_t *arg );
void hqr_high_greedy_init   ( hqr_subpiv_t *arg, int minMN );
void hqr_high_greedy1p_init ( hqr_subpiv_t *arg );
/**
 * @}
 *
 * @name Matrix reduction trees
 * @{
 */

void libhqr_fct_to_mtx( const libhqr_tree_t *in, libhqr_tree_t *out );
void libhqr_rdmtx_initfct( libhqr_tree_t *out );

/**
 * @}
 */

#endif /* _libhqr_internal_h_ */
