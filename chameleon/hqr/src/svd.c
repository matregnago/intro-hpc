/**
 *
 * @file svd.c
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 * This file contains all the function to describe the dependencies
 * used in the Xgeqrf_param.jdf file.
 * The QR factorization done with this file relies on three levels:
 *     - the first one is using a flat tree with TS kernels. The
 *       height of this tree is defined by the parameter 'a'. If 'a'
 *       is set to A->mt, the factorization is identical to the one
 *       perform by PLASMA_zgeqrf.
 *       For all subdiagonal "macro-tiles", the line reduced is always
 *       the first.  For all diagonal "macro-tiles", the factorization
 *       performed is identical to the one performed by PLASMA_zgeqrf.
 *
 *     - the third level is using a reduction tree of size 'p'. By
 *       default, the parameter 'p' should be equal to the number of
 *       processors used for the computation, but can be set
 *       differently. (see further example). The type of tree used at
 *       this level is defined by the hlvl parameter. It can be flat
 *       or greedy.
 *       CODE DETAILS: This tree and all the function related to it
 *       are performing a QR factorization on a band matrix with 'p'
 *       the size of the band. All the functions take global indices
 *       as input and return global indices as output.
 *
 *     - Finally, a second 'low' level of reduction tree is applied.
 *       The size of this tree is induced by the parameters 'a' and 'p'
 *       from the first and third levels and is A->mt / ( p * a ). This
 *       tree is reproduced p times for each subset of tiles
 *       S_k = {i in [0, A->mt-1] \ i%p*a = k } with k in [0, p-1].
 *       The tree used for the reduction is defined by the llvl
 *       parameter and can be: flat, greedy, fibonacci or binary.
 *       CODE DETAILS: For commodity, the size of this tree is always
 *       ceil(A->mt / (p * a) ) inducing some extra tests in the code.
 *       All the functions related to this level of tree take as input
 *       the local indices in the A->mt / (p*a) matrix and the global
 *       k. They return the local index. The reductions are so
 *       performed on a trapezoidal matrices where the step is defined
 *       by a:
 *                                    <- min( lhlvl_mt, min( mt, nt ) ) ->
 *                                     __a__   a     a
 *                                    |     |_____
 *                                    |           |_____
 *                                    |                 |_____
 *        llvl_mt = ceil(MT/ (a*p))   |                       |_____
 *                                    |                             |_____
 *                                    |___________________________________|
 *
 *
 *
 *   At each step of the factorization, the lines of tiles are divided
 *   in 4 types:
 *     - QRPARAM_TILE_TS: They are the lines annihilated by a TS
 *     kernel, these lines are never used as an annihilator.  They are
 *     the lines i, with 1 < (i/p)%a < a and i > (k+1)*p
 *     - QRPARAM_TILE_LOCALTT: They are the lines used as annhilitor
 *     in the TS kernels annihiling the QRPARAM_TILE_TS lines.  They
 *     are themselves annihilated by the TT kernel of the low level
 *     reduction tree.  The roots of the local trees are the lines i,
 *     with i/p = k.
 *     - QRPARAM_TILE_DOMINO: These are the lines that are
 *     annhilihated with a domino effect in the band defined by (i/p)
 *     <= k and i >= k
 *     - QRPARAM_TILE_DISTTT: These are the lines annihilated by the
 *     high level tree to reduce communications.
 *     These lines are defined by (i-k)/p = 0.
 */
#include "libhqr_internal.h"
#include <stdlib.h>
#include <string.h>

/*
 * Common functions
 */
static int svd_getnbgeqrf( const libhqr_tree_t *qrtree, int k );
static int svd_getm(       const libhqr_tree_t *qrtree, int k, int i   );
static int svd_geti(       const libhqr_tree_t *qrtree, int k, int m   );
static int svd_gettype(    const libhqr_tree_t *qrtree, int k, int m   );

#define svd_getipiv( __qrtree, _k ) ((__qrtree)->llvl->ipiv + ((__qrtree)->llvl->ldd) * (_k) )
#define svd_geta( __qrtree, _k ) ( (svd_getipiv( (__qrtree), (_k) ))[0] )

/*
 * Extra parameter:
 *    gmt - Global number of tiles in a column of the complete distributed matrix
 * Return:
 *    The number of geqrt to execute in the panel k
 */
static int
svd_getnbgeqrf( const libhqr_tree_t *qrtree,
                int k )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int p = qrtree->p;
    int gmt = qrtree->mt;
    int a = svd_geta(arg, k);
    int pa = p * a;
    int nb_1, nb_2, nb_3;
    int nb_11, nb_12;

    /* Number of tasks of type 3 */
    nb_3 = p;

    /* Number of extra tile of type 1 between the tile of type 3 and the first of nb11 */
    nb_2 = nbextra1_formula;

    /* First multiple of p*a under the diagonal of step 1 */
    nb_11 = ( (k + p + pa - 1 ) / pa ) * pa;

    /* Last multiple of p*a lower than A->mt */
    nb_12 = ( gmt / pa ) * pa;

    /* Number of tasks of type 1 between nb_11 and nb_12 */
    nb_1 = (nb_12 - nb_11) / a;

    /* Add leftover */
    nb_1 += libhqr_imin( p, gmt - nb_12 );

    return libhqr_imin( nb_1 + nb_2 + nb_3, gmt - k);
}

/*
 * Extra parameter:
 *    i - indice of the geqrt in the continuous space
 * Return:
 *    The global indice m of the i th geqrt in the panel k
 */
static int
svd_getm( const libhqr_tree_t *qrtree,
          int k, int i )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int p = qrtree->p;
    int a = svd_geta(arg, k);

    int pa = p * a;
    int nbextra1 = nbextra1_formula;
    int nb23 = p + nbextra1;

    /* Tile of type 2 or 3 or the 1 between the diagonal and the multiple after the diagonal */
    if ( i < nb23 ) {
        return k+i;
    }
    /* Tile of type 1 */
    else {
        int j = i - nb23;
        int pos1 = ( ( (p + k    ) + pa - 1 ) / pa ) * pa;
        return pos1 + (j/p) * pa + j%p;
    }
}

/*
 * Extra parameter:
 *    m - The global indice m of a geqrt in the panel k
 * Return:
 *    The index i of the geqrt in the panel k
 */
static int
svd_geti( const libhqr_tree_t *qrtree,
          int k, int m )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int p = qrtree->p;
    int a = svd_geta(arg, k);

    int pa = p * a;
    int nbextra1 = nbextra1_formula;
    int end2 = p + k + nbextra1;

    /* Tile of type 2 or 3 or the 1 between the diagonal and the multiple after the diagonal */
    if ( m < end2 ) {
        return m-k;
    }
    /* Tile of type 1 */
    else {
        int pos1 = ( ( (p + k) + pa - 1 ) / pa ) * pa;
        int j    = m - pos1;
        int nb23 = p + nbextra1;
        return nb23 + (j / pa) * p + j%pa;
    }
}

/*
 * Extra parameter:
 *      m - The global indice m of the row in the panel k
 * Return
 *     -1 - Error
 *      0 - if m is reduced thanks to a TS kernel
 *      1 - if m is reduced thanks to the low level tree
 *      2 - if m is reduced thanks to the bubble tree
 *      3 - if m is reduced in distributed
 */
static int
svd_gettype( const libhqr_tree_t *qrtree,
             int k, int m )
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int p = qrtree->p;
    int a = svd_geta(arg, k);

    /* Element to be reduce in distributed */
    if (m < k + p) {
        return 3;
    }
    /* Lower triangle of the matrix */
    else {
        if( (m / p) % a == 0 ) {
            return 1;
        }
        else {
            return 0;
        }
    }
}

/****************************************************
 *
 *   Generic functions currpiv,prevpiv,nextpiv
 *
 ***************************************************/
static int
svd_currpiv(const libhqr_tree_t *qrtree, int k, int m)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, tmpk;
    int lm, rank;
    int a = svd_geta( arg, k );
    int p = qrtree->p;
    int gmt = qrtree->mt;

    lm   = m / p; /* Local index in the distribution over p domains */
    rank = m % p; /* Staring index in this distribution             */

    /* TS level common to every case */
    switch( svd_gettype( qrtree, k, m ) )
    {
    case 0:
        tmp = lm / a;
        /* tmpk = (k + p - 1 - m%p) / p / a;  */
        tmpk = k / (p * a);
        return ( tmp == tmpk ) ? k + (m-k)%p : tmp * a * p + rank;
        break;
    case 1:
        tmp = arg->llvl->currpiv(arg->llvl, k, m);
        /* tmpk = (k + p - 1 - m%p) / p / a; */
        tmpk = k / (p * a);
        return ( tmp == tmpk ) ? k + (m-k)%p : tmp * a * p + rank;
        break;
    case 2:
        assert(0);
        break;
    case 3:
        if ( arg->hlvl != NULL ) {
            return arg->hlvl->currpiv(arg->hlvl, k, m);
        }
        libhqr_attr_fallthrough;
    default:
        return gmt;
    }
    return -1;
}

/**
 *  svd_nextpiv - Computes the next row killed by the row p, after
 *  it has kill the row start.
 *
 * @param[in] k
 *         Factorization step
 *
 * @param[in] pivot
 *         Line used as killer
 *
 * @param[in] start
 *         Starting point to search the next line killed by p after start
 *         start must be equal to A.mt to find the first row killed by p.
 *         if start != A.mt, start must be killed by p
 *
 * @return:
 *   - -1 if start doesn't respect the previous conditions
 *   -  m, the following row killed by p if it exists, A->mt otherwise
 */
static int svd_nextpiv(const libhqr_tree_t *qrtree, int k, int pivot, int start)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, ls, lp, nextp;
    int rpivot, lstart;
    int *ipiv = svd_getipiv( arg, k );
    int a = ipiv[0];
    int ldd = ipiv[1];
    int p = qrtree->p;
    int gmt = qrtree->mt;

    rpivot = pivot % p; /* Staring index in this distribution             */

    /* Local index in the distribution over p domains */
    lstart = ( start == gmt ) ? ldd * a : start / p;

    myassert( start > pivot && pivot >= k );
    myassert( start == gmt || pivot == svd_currpiv( qrtree, k, start ) );

    /* TS level common to every case */
    ls = (start < gmt) ? svd_gettype( qrtree, k, start ) : -1;
    lp = svd_gettype( qrtree, k, pivot );

    switch( ls )
    {
    case LIBHQR_KILLED_BY_DOMINO:
        assert(0);
        libhqr_attr_fallthrough;

    case -1:

        if ( lp == LIBHQR_KILLED_BY_TS ) {
            myassert( start == gmt );
            return gmt;
        }
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_TS:
        if ( start == gmt ) {
            nextp = pivot + p;
        }
        else {
            nextp = start + p;
        }

        if ( ( nextp < gmt ) &&
             ( nextp < (pivot + a*p) ) &&
             ( (nextp/p)%a != 0 ) )
        {
            return nextp;
        }
        //start = gmt;  /* Commented to silent SonarQube */
        lstart = ldd * a;
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_LOCALTREE:
        /* Get the next pivot for the low level tree */
        tmp = arg->llvl->nextpiv(arg->llvl, k, pivot, lstart / a );

        if ( ((tmp * a * p + rpivot) >= gmt) &&
             (tmp == ldd-1) )
        {
            tmp = arg->llvl->nextpiv(arg->llvl, k, pivot, tmp);
        }

        if ( tmp != ldd ) {
            return tmp * a * p + rpivot;
        }

        /* no next of type 1, we reset start to search the next 2 */
        start = gmt;
        /* lstart = ldd * a; */
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_DISTTREE:

        if ( lp < LIBHQR_KILLED_BY_DISTTREE ) {
            return gmt;
        }

        if( arg->hlvl != NULL ) {
            tmp = arg->hlvl->nextpiv( arg->hlvl, k, pivot, start );
            if ( tmp != gmt ) {
                return tmp;
            }
        }
        libhqr_attr_fallthrough;

    default:
        return gmt;
    }
}

/**
 *  svd_prevpiv - Computes the previous row killed by the row p, before
 *  to kill the row start.
 *
 * @param[in] k
 *         Factorization step
 *
 * @param[in] pivot
 *         Line used as killer
 *
 * @param[in] start
 *         Starting point to search the previous line killed by p before start
 *         start must be killed by p, and start must be greater or equal to p
 *
 * @return:
 *   - -1 if start doesn't respect the previous conditions
 *   -  m, the previous row killed by p if it exists, A->mt otherwise
 */
static int
svd_prevpiv(const libhqr_tree_t *qrtree, int k, int pivot, int start)
{
    hqr_args_t *arg = (hqr_args_t*)(qrtree->args);
    int tmp, ls, lp, nextp;
    int lpivot, rpivot, lstart;
    int *ipiv = svd_getipiv( arg, k );
    int a = ipiv[0];
    int ldd = ipiv[1];
    int p = qrtree->p;
    int gmt = qrtree->mt;

    lpivot = pivot / p; /* Local index in the distribution over p domains */
    rpivot = pivot % p; /* Staring index in this distribution             */
    lstart = start / p; /* Local index in the distribution over p domains */

    myassert( start >= pivot && pivot >= k && start < gmt );
    myassert( start == pivot || pivot == svd_currpiv( qrtree, k, start ) );

    /* TS level common to every case */
    ls = svd_gettype( qrtree, k, start );
    lp = svd_gettype( qrtree, k, pivot );

    if ( lp == LIBHQR_KILLED_BY_TS ) {
        return gmt;
    }

    myassert( lp >= ls );
    switch( ls )
    {
    case LIBHQR_KILLED_BY_DOMINO:
        assert(0);
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_DISTTREE:
        if( arg->hlvl != NULL ) {
            tmp = arg->hlvl->prevpiv( arg->hlvl, k, pivot, start );
            if ( tmp != gmt ) {
                return tmp;
            }
        }

        //start = pivot; /* Commented to silent SonarQube */
        lstart = pivot / p;
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_LOCALTREE:
        tmp = arg->llvl->prevpiv(arg->llvl, k, pivot, lstart / a);

        if ( ((tmp * a * p + rpivot) >= gmt) &&
             (tmp == ldd-1) )
        {
            tmp = arg->llvl->prevpiv(arg->llvl, k, pivot, tmp);
        }

        if ( tmp != ldd ) {
            return tmp * a * p + rpivot;
        }

        start = pivot;
        libhqr_attr_fallthrough;

    case LIBHQR_KILLED_BY_TS:
        if ( start == pivot ) {
            tmp = lpivot + a - 1 - lpivot%a;
            nextp = tmp * p + rpivot;

            while( (pivot < nextp) &&
                   (nextp >= gmt ) )
            {
                nextp -= p;
            }
        } else {
            nextp = start - p; /*(lstart - 1) * p + rpivot;*/
        }
        assert(nextp < gmt);
        if ( pivot < nextp ) {
            return nextp;
        }
        libhqr_attr_fallthrough;

    default:
        return gmt;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup libhqr
 *
 * libhqr_svd_init - Create the tree structures that will describes the
 * operation performed during QR/LQ reduction step of the gebrd_ge2gb operation.
 *
 * Trees available parameters are described below. It is recommended to:
 *   - set p to the same value than the P-by-Q process grid used to distribute
 *     the data. (P for QR factorization, Q for LQ factorization).
 *   - set the low level tree to LIBHQR_GREEDY_TREE.
 *   - set the high level tree to:
 *         1) LIBHQR_FLAT_TREE when the problem is square, because it divides
 *            by two the volume of communication of any other tree.
 *         2) LIBHQR_FIBONACCI_TREE when the problem is tall and skinny (QR) or
 *            small and fat (LQ), because it reduces the critical path length.
 *   - Disable the domino effect when problem is square, to keep high efficiency
 *     kernel proportion high.
 *   - Enable the domino effect when problem is tall and skinny (QR) or
 *     small and fat (LQ) to increase parallelism and reduce critical path length.
 *   - Round-robin on TS domain (tsrr) option should be disabled. It is
 *     experimental and is not safe.
 *
 * These are the default when a parameter is set to -1;
 *
 * See http://www.netlib.org/lapack/lawnspdf/lawn257.pdf
 *
 *******************************************************************************
 *
 * @param[in,out] qrtree
 *          On entry, an allocated structure uninitialized.
 *          On exit, the structure initialized according to the given parameters.
 *
 * @param[in] trans
 *          @arg LIBHQR_QR:   Structure is initialized for the QR steps.
 *          @arg LIBHQR_LQ:     Structure is initialized for the LQ steps.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A to be factorized, on which
 *          QR/LQ reduction steps will be performed. In case, of
 *          R-bidiagonalization, don't forget to provide the square submatrix
 *          that is concerned by those operations.
 *          The descriptor is untouched and only mt/nt/P parameters are used.
 *
 * @param[in] type_hlvl
 *          Defines the tree used to reduce the main tiles of each domain. This
 *          is a band lower diagonal matrix of width p.
 *          @arg LIBHQR_FLAT_TREE: A Flat tree is used to reduce the tiles.
 *          @arg LIBHQR_GREEDY_TREE: A Greedy tree is used to reduce the tiles.
 *          @arg LIBHQR_FIBONACCI_TREE: A Fibonacci tree is used to reduce the
 *          tiles.
 *          @arg LIBHQR_BINARY_TREE: A Binary tree is used to reduce the tiles.
 *          @arg LIBHQR_GREEDY1P_TREE: A Greedy tree is computed for the first
 *          column and then duplicated on all others.
 *          @arg -1: The default is used (LIBHQR_FIBONACCI_TREE)
 *
 * @param[in] p
 *          Defines the number of distributed domains, ie the width of the high
 *          level reduction tree.  If p == 1, no high level reduction tree is
 *          used. If p == mt, this enforce the high level reduction tree to be
 *          performed on the full matrix.
 *          By default, it is recommended to set p to P if trans ==
 *          LIBHQR_QR, and to Q otherwise, where P-by-Q is the process grid
 *          used to distributed the data. (p > 0)
 *
 * @param[in] nbthread_per_node
 *          Define the number of working threads per node to configure the
 *          adaptativ local tree to provide at least (ratio * nbthread_per_node)
 *          tasks per step when possible by creating the right amount of TS and
 *          TT kernels.
 *
 * @param[in] ratio
 *          Define the minimal number of tasks per thread that the adaptiv tree
 *          must provide at the lowest level of the tree.
 *
 *******************************************************************************
 *
 * @retval -i if the ith parameters is incorrect.
 * @retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa libhqr_hqr_finalize
 * @sa libhqr_hqr_init
 * @sa dplasma_zgeqrf_param
 * @sa dplasma_cgeqrf_param
 * @sa dplasma_dgeqrf_param
 * @sa dplasma_sgeqrf_param
 *
 ******************************************************************************/
int
libhqr_initfct_svd( libhqr_tree_t *qrtree,
                    libhqr_facto_e trans, const libhqr_matrix_t *A,
                    int type_hlvl, int p, int nbthread_per_node, int ratio )
{
    int low_mt, minMN, a = -1;
    hqr_args_t *arg;

    if (qrtree == NULL) {
        fprintf(stderr,"libhqr_initfct_svd: illegal value of qrtree");
        return -1;
    }
    if ((trans != LIBHQR_QR) &&
        (trans != LIBHQR_LQ)) {
        fprintf(stderr, "libhqr_initfct_svd: illegal value of trans");
        return -2;
    }
    if (A == NULL) {
        fprintf(stderr, "libhqr_initfct_svd: illegal value of A");
        return -3;
    }

    /* Compute parameters */
    p = libhqr_imax( p, 1 );
    minMN = libhqr_imin(A->mt, A->nt);

    memset( qrtree, 0, sizeof( libhqr_tree_t ) );
    qrtree->getnbgeqrf = svd_getnbgeqrf;
    qrtree->getm       = svd_getm;
    qrtree->geti       = svd_geti;
    qrtree->gettype    = svd_gettype;
    qrtree->currpiv    = svd_currpiv;
    qrtree->nextpiv    = svd_nextpiv;
    qrtree->prevpiv    = svd_prevpiv;

    qrtree->init  = LIBHQR_QRTREE_SVD;
    qrtree->facto = trans;
    qrtree->mt    = (trans == LIBHQR_QR) ? A->mt : A->nt;
    qrtree->nt    = minMN;

    qrtree->a    = a;
    qrtree->p    = p;
    qrtree->args = NULL;

    arg = (hqr_args_t*) malloc( sizeof(hqr_args_t) );

    arg->llvl = (hqr_subpiv_t*) malloc( sizeof(hqr_subpiv_t) );
    arg->hlvl = NULL;

    low_mt = (qrtree->mt + p - 1) / ( p );

    arg->llvl->minMN  = minMN;
    arg->llvl->ldd    = low_mt;
    arg->llvl->a      = a;
    arg->llvl->p      = p;
    arg->llvl->domino = 0;

    svd_low_adaptiv_init( arg->llvl, qrtree->mt, qrtree->nt,
                          nbthread_per_node * (A->nodes / p), ratio );

    if ( p > 1 ) {
        arg->hlvl = (hqr_subpiv_t*) malloc( sizeof(hqr_subpiv_t) );

        arg->llvl->minMN  = minMN;
        arg->hlvl->ldd    = qrtree->mt;
        arg->hlvl->a      = a;
        arg->hlvl->p      = p;
        arg->hlvl->domino = 0;

        switch( type_hlvl ) {
        case LIBHQR_FLAT_TREE :
            hqr_high_flat_init(arg->hlvl);
            break;
        case LIBHQR_GREEDY_TREE :
            hqr_high_greedy_init(arg->hlvl, minMN);
            break;
        case LIBHQR_GREEDY1P_TREE :
            hqr_high_greedy1p_init(arg->hlvl);
            break;
        case LIBHQR_BINARY_TREE :
            hqr_high_binary_init(arg->hlvl);
            break;
        case LIBHQR_FIBONACCI_TREE :
            hqr_high_fibonacci_init(arg->hlvl);
            break;
        default:
            hqr_high_fibonacci_init(arg->hlvl);
        }
    }

    qrtree->args = (void*)arg;

    return 0;
}

int
libhqr_initmtx_svd( libhqr_tree_t *qrtree,
                    libhqr_facto_e trans, const libhqr_matrix_t *A,
                    int type_hlvl, int p, int nbthread_per_node, int ratio )
{
    libhqr_tree_t qrtreefct;
    int rc;

    /* Missing functions for now, should not be used */
    assert(0);

    rc = libhqr_initfct_svd( &qrtreefct, trans, A, type_hlvl, p, nbthread_per_node, ratio );
    if ( rc != 0 ) {
        return rc;
    }

    libhqr_fct_to_mtx( &qrtreefct, qrtree );

    /* Free the initial qrtree */
    libhqr_finalize( &qrtreefct );

    return 0;
}

int
libhqr_init_svd( libhqr_tree_t *qrtree,
                 libhqr_facto_e trans, const libhqr_matrix_t *A,
                 int type_hlvl, int p, int nbthread_per_node, int ratio )
{
    return libhqr_initfct_svd( qrtree, trans, A, type_hlvl, p, nbthread_per_node, ratio );
}
