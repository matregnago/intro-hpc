/**
 *
 * @file core_ipiv_to_perm.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_ipiv_to_perm CPU kernel
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @date 2025-12-19
 */
#include "coreblas.h"

/**
 *******************************************************************************
 *
 * The idea here is to generate a permutation from the sequence of
 * pivot.  To avoid storing one whole column at each step, we keep
 * track of two vectors of nb elements, the first one contains the
 * permutation of the first nb elements, and the second one contains
 * the inverse permutation of those same elements.
 *
 * Lets have i the element to pivot with ip. ipiv[i] = ip;
 * We set i_1 as such invp[ i_1  ] = i
 *  and  ip_1 as such invp[ ip_1 ] = ip
 *
 * At each step we want to:
 *   - swap perm[i] and perm[ip]
 *   - set invp[i_1] to ip
 *   - set invp[ip_1] to i
 *
 *******************************************************************************
 *
 * @param[in] m0
 *          The base index for all values in ipiv, perm and invp. m0 >= 0.
 *
 * @param[in] m
 *          The number of elements in perm and invp. m >= 0.
 *
 * @param[in] k
 *          The number of elements in ipiv. k >= 0.
 *
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done.
 *
 * @param[in] ipiv
 *          The pivot array of size n. This is a (m0+1)-based indices array to follow
 *          the Fortran standard.
 *
 * @param[out] perm
 *          The permutation array of the destination row indices (m0-based) of the [1,n] set of rows.
 *
 * @param[out] invp
 *          The permutation array of the origin row indices (m0-based) of the [1,n] set of rows.
 *
 */
void CORE_ipiv_to_perm( int m0, int m, int k, int K1, int K2,
                        const int *ipiv, int *perm, int *invp )
{
    /*
     * Let's note that variables suffixed by l are local values (shifted by
     * -m0), and variable suffixed by g are global values.
     */
    int il, ipl, il_1, ipl_1, jl;
    int ig, ipg, ig_1, ipg_1;

    /* Loop through perm and invp to initialize them with no pivoting */
    for( il=0; il < m; il++ ) {
        ig = il + m0;
        perm[il] = ig;
        invp[il] = ig;
    }

    /* Loop through ipiv to compute perm and invp */
    for(il = 0; il < k; il++) {
        ig = il + m0;
        if ( ( ig < K1 ) || ( ig > K2 ) ) {
            continue;
        }

        ipg = ipiv[ il ] - 1;
        ipl = ipg - m0;

        /* Pivot should only returns rows below or equal to the current one */
        assert( ipl >= il );

        /* If the row i is permuted with the row ip (i != ip) */
        if ( ipl > il ) {

            /* Save initial permutation in the cycle */
            ig_1 = perm[ il ];

            /* If the row ipl is in the current block */
            if ( ipl < m ) {
                /* Save final permutation point in the cycle and update intermediate one */
                ipg_1 = perm[ ipl ];
                perm[ ipl ] = ig_1;
            }
            /* If the row ipl is not in the current block */
            else {
                ipg_1 = ipg;
                /* Loop through invp to see if ip has already been set in the current block */
                for( jl=0; jl < m; jl++ ) {
                    if( invp[jl] == ipg ) {
                        ipg_1 = jl + m0;
                        break;
                    }
                }
            }

            /* Save the final permutation */
            perm[ il ] = ipg_1;
            il_1  = ig_1  - m0;
            ipl_1 = ipg_1 - m0;

            /* Update initial and final invp if necessary */
            if (il_1  < m) invp[il_1 ] = ipg;
            if (ipl_1 < m) invp[ipl_1] = ig;
        }
    }
}

/**
 * @brief Convert a perm array to a perm_idx one. See CORE_perm_to_idx for
 * descriptions.
 */
static inline void
convert_one_perm_array( int        m0,
                        int        m,
                        int        mt,
                        const int *perm,
                        int       *perm_idx )
{
    int  i, ipl, ipm, ipg, idx;
    int  cnt_p[ 2 * mt ];
    int *perm_data = perm_idx + mt + 1;
    memset( cnt_p, 0, sizeof(int) * mt );

    /* Intialize the first tile */
    cnt_p[ 0 ] = m;

    /* Loop through perm and invp to count the number of lines to get/set per tile */
    for( i = 0; i < m; i++ ) {
        /* Update number of lines to get for perm */
        ipl = perm[ i ] - m0;
        ipm = ipl / m;

        if ( ipm > 0 ) {
            cnt_p[ 0   ]--;
            cnt_p[ ipm ]++;
        }
    }

    /* Udpate index array */
    perm_idx[ 0 ] = 0;
    for( i=0; i < mt; i++ ) {
        perm_idx[ i+1 ] = perm_idx[ i ] + cnt_p[ i ];
        cnt_p[ i ] = 0;
    }

    /* Store the couples sorted by tile */
    for( i = 0; i < m; i++ ) {
        ipg = perm[i];
        ipl = ipg - m0;
        ipm = ipl / m;
        idx = perm_idx[ ipm ] + cnt_p[ ipm ];
        cnt_p[ ipm ]++;

        perm_data[ idx * 2 ]     = i;
        perm_data[ idx * 2 + 1 ] = ipg;
    }
}

/**
 *******************************************************************************
 *
 * This function converts the perm and invp array to a structure similar to a
 * "csr" that list the rows and their values sorted by tile.
 *
 * For a given permutation, the perm_idx array will consists of two arrays, the
 * first one, similar to colptr in csr, that gives the indices of the elements
 * for each tile.
 * The second one is an array of couples (local row index in the final tile
 * after permutation, global row index of the selected row)
 *
 * Example:
 * idxptr = { 0, 5, 9, 13 }
 * data   = {  2  3  4  5  9  0  7  8 12  1  6 10 11
 *            20 23 25 21 16 26 32 38 30 48 41 40 43 }
 * perm_idx stores idxptr immediately followed by data in column-major.
 *
 *******************************************************************************
 *
 * @param[in] m0
 *          The base index for all values in ipiv, perm and invp. m0 >= 0.
 *
 * @param[in] m
 *          The number of elements in perm and invp. m >= 0.
 *
 * @param[in] mt
 *          The number of tiles that stores all those possible lines. mt >= 1.
 *
 * @param[in] perm
 *          The permutation array generated by CORE_ipiv_to_perm() of size m.
 *          It contains the destination row indices (0-based) of the [1,m] set of rows
 *
 * @param[in] invp
 *          The inverse permutation array generated by CORE_ipiv_to_perm() of size m.
 *          It contains the origin row indices (0-based) of the [1,m] set of rows.
 *
 * @param[out] perm_idx
 *          The permutation array index with the number of rows per tile and the
 *          permutations couples sorted by tile.
 *
 * @param[out] invp_idx
 *          The inverse permutation array index with the number of rows per tile and the
 *          permutations couples sorted by tile.
 *
 */
void CORE_perm_to_idx( int        m0,
                       int        m,
                       int        mt,
                       const int *perm,
                       const int *invp,
                       int       *perm_idx,
                       int       *invp_idx )
{
    convert_one_perm_array( m0, m, mt, perm, perm_idx );
    convert_one_perm_array( m0, m, mt, invp, invp_idx );
}
