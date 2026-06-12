/**
 *
 * @file descriptor_helpers.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @author Lionel Eyraud-Dubois
 * @author Alycia Lisito
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @date 2025-12-19
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users
 *
 */
#include "control/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "control/descriptor.h"
#include "chameleon/runtime.h"

/**
 * @brief Return the rank of the tile A( m, n ) in a classic 2D Block Cyclic
 * distribution PxQ.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The rank of the tile A( m, n )
 *
 */
int chameleon_getrankof_2d( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    int nn = n + A->j / A->nb;
    return (mm % chameleon_desc_datadist_get_iparam(A,0)) * chameleon_desc_datadist_get_iparam(A,1) + (nn % chameleon_desc_datadist_get_iparam(A,1));
}

/**
 * @brief Return the rank associated to the diagonal tile ( m, m ) of a classic
 * 2D Block Cyclic distribution PxQ.
 *
 * @param[in] A
 *        A specific matrix descriptor that holds only diagonal tiles. Thus, n is never used.
 *
 * @param[in] m
 *        The row and column index of the tile.
 *
 * @param[in] n
 *        Unused
 *
 * @return The rank of the tile A( m, m )
 *
 */
int chameleon_getrankof_2d_diag( const CHAM_desc_t *A, int m, int n )
{
    int mm = m + A->i / A->mb;
    (void)n;
    return (mm % chameleon_desc_datadist_get_iparam(A,0)) * chameleon_desc_datadist_get_iparam(A,1) + (mm % chameleon_desc_datadist_get_iparam(A,1));
}

/**
 * @brief Return the rank of the process responsible for the permutation of the tile (m, n)
 * in a classic 2D Block Cyclic distribution PxQ.
 *
 * @param[in] IPIV
 *        The ipiv descriptor.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The rank of the process responsible for the row permutation of the tile (m, n)
 *
 */
int chameleon_getrankof_ipiv_2d_row( const CHAM_ipiv_t *IPIV, int m, int n )
{
    int Q = IPIV->NP / IPIV->P;
    (void)n;
    return ( m % IPIV->P ) * Q;
}

/**
 * @brief Return the rank of the process responsible for the column permutation of the tile (m, n)
 * in a classic 2D Block Cyclic distribution PxQ.
 *
 * @param[in] IPIV
 *        The ipiv descriptor.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The rank of the process responsible for the permutation of the tile (m, n)
 *
 */
int chameleon_getrankof_ipiv_2d_col( const CHAM_ipiv_t *IPIV, int m, int n )
{
    int Q = IPIV->NP / IPIV->P;
    (void)m;
    return n % Q;
}

/**
 * @brief Return the rank of the process responsible for the permutation of the tile (m, n)
 * when used for getrf in a classic 2D Block Cyclic distribution PxQ.
 *
 * @param[in] IPIV
 *        The ipiv descriptor.
 *
 * @param[in] m
 *        The row and column index of the tile.
 *
 * @param[in] n
 *        Unused
 *
 * @return The rank of the process responsible for the permutation of the tile (m, n)
 *
 */
int chameleon_getrankof_ipiv_2d_diag( const CHAM_ipiv_t *IPIV, int m, int n )
{
    (void)n;
    int Q = IPIV->NP / IPIV->P;
    return (m % IPIV->P) * Q + (m % Q);
}

/**
 * @brief Test if the current MPI process is involved in the panel k for 2DBC distributions.
 *
 * @param[in] A
 *        The matrix descriptor.
 *
 * @param[in] k
 *        The index of the panel to test.
 *
 * @return 1 if the current MPI process contributes to the panel k.
 *         0 if the current MPI process doesn't contribute to the panel k.
 *
 */
int chameleon_involved_in_panelk_2dbc( const CHAM_desc_t *A, int k ) {
    int myrank = A->myrank;
    return ( myrank % chameleon_desc_datadist_get_iparam(A,1) == k % chameleon_desc_datadist_get_iparam(A,1) );
}

/**
 * @brief Test if the MPI process p is involved in the panel k for 2DBC distributions.
 *
 * @param[in] A
 *        The matrix descriptor.
 *
 * @param[in] k
 *        The index of the panel to test.
 *
 * @param[in] p
 *        The rank of the MPI process.
 *
 * @return 1 if the current MPI process contributes to the panel k.
 *         0 if the current MPI process doesn't contribute to the panel k.
 *
 */
int chameleon_p_involved_in_panelk_2dbc( const CHAM_desc_t *A, int k, int p ) {
    return ( p % chameleon_desc_datadist_get_iparam(A,1) == k % chameleon_desc_datadist_get_iparam(A,1) );
}

/**
 * @brief Test if the current MPI process is involved in the panel k for 2DBC distributions.
 *
 * @param[in] A
 *        The matrix descriptor.
 *
 * @param[in] k
 *        The index of the panel to test.
 *
 * @param[in] n
 *        The index of the panel to test.
 *
 * @param[inout] ws_reduce
 *        The i.
 *
 */
void chameleon_get_proc_involved_in_panelk_2dbc( const CHAM_desc_t *A,
                                                 int                k,
                                                 int                n,
                                                 void              *ws_reduce )
{
#if defined (CHAMELEON_USE_MPI)
    CHAM_reduce_t *reduce        = (CHAM_reduce_t*) ws_reduce;
    int           *proc_involved = reduce->proc_involved;
    int            b, rank, np;

    np = 0;
    reduce->involved = 0;
    for ( b = k; (b < A->mt) && ((b-k) < chameleon_desc_datadist_get_iparam(A, 0)); b ++ ) {
        rank = A->get_rankof( A, b, n );
        proc_involved[ b-k ] = rank;
        np ++;
        if ( rank == A->myrank ) {
            reduce->involved = 1;
        }
    }
    reduce->proc_involved = proc_involved;
    reduce->np_involved   = np;
#else
    (void)A;
    (void)k;
    (void)n;
    (void)ws_reduce;
#endif
}

/**
 * @brief Test if the current MPI process is involved in the panel k for 2DBC distributions.
 *
 * @param[in] A
 *        The matrix descriptor.
 *
 * @param[in] m
 *        The index of the panel to test.
 *
 * @param[in] k
 *        The index of the panel to test.
 *
 * @param[inout] ws_reduce
 *        The i.
 *
 */
void chameleon_get_proc_involved_in_rowpanelk_2dbc( const CHAM_desc_t *A,
                                                    int                m,
                                                    int                k,
                                                    void              *ws_reduce )
{
#if defined (CHAMELEON_USE_MPI)
    CHAM_reduce_t *reduce = (CHAM_reduce_t*) ws_reduce;
    int           *proc_involved = reduce->proc_involved;
    int            b, rank, np;

    np = 0;
    reduce->involved = 0;
    for ( b = k; (b < A->nt) && ((b-k) < chameleon_desc_datadist_get_iparam(A, 1)); b ++ ) {
        rank = chameleon_getrankof_2d( A, m, b );
        proc_involved[ b-k ] = rank;
        np ++;
        if ( rank == A->myrank ) {
            reduce->involved = 1;
        }
    }
    reduce->proc_involved = proc_involved;
    reduce->np_involved   = np;
#else
    (void)A;
    (void)k;
    (void)m;
    (void)ws_reduce;
#endif
}

/**
 * @brief Initializes a custom distribution based on an external file.
 *
 *  External file format: First line contains the dimensions M and N (space-separated)
 *  Next M lines each have N integers (with values from 1 to number of nodes)
 *
 *
 * @param[out] custom_dist
 *        On exit, the resulting custom distribution
 *
 * @param[in] dist_file
 *        The path to the external file to be read
 *
 * @return CHAMELEON_SUCCESS on successful exit, CHAMELEON_ERR_OUT_OF_RESOURCES
 *         or CHAMELEON_ERR_ILLEGAL_VALUE on issue.
 *
 */
int chameleon_getrankof_custom_init( custom_dist_t **custom_dist,
                                     const char     *dist_file )
{
    custom_dist_t *result;
    FILE          *f;
    int            i, j, rc;
    int            np, dist_m, dist_n;

    *custom_dist = NULL;

    /* Get number of processes to check for correctness */
    np = CHAMELEON_Comm_size();

    /* Allocate memory */
    result = (custom_dist_t*) malloc( sizeof(custom_dist_t) );
    if ( result == NULL ) {
        chameleon_error( "chameleon_getrankof_custom_init", "malloc() failed" );
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    result->dist_file = dist_file;
    f = fopen( result->dist_file , "r" );
    if ( f == NULL ) {
        char message[300];
        snprintf( message, 300, "could not open file '%s'", dist_file );
        chameleon_error( "chameleon_getrankof_custom_init", message );
        free( result );
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    rc = fscanf( f, "%d %d", &dist_m, &dist_n );
    if ( rc < 2 ) {
        char message[300];
        snprintf( message, 300, "could not read m and n in file '%s'", dist_file );
        chameleon_error( "chameleon_getrankof_custom_init", message );
        free( result );
        fclose( f );
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    if ( (dist_m <= 0 ) || (dist_n <= 0) ) {
        char message[300];
        snprintf( message, 300, "Incorrect values for dist_m(%d) and/or dist_n(%d)", dist_m, dist_n );
        chameleon_error( "chameleon_getrankof_custom_init", message );
        free( result );
        fclose( f );
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    result->dist_m = dist_m;
    result->dist_n = dist_n;

    result->blocks_dist = (int*) malloc( sizeof(int) * dist_m * dist_n );
    if ( result->blocks_dist == NULL ) {
        chameleon_error( "chameleon_getrankof_custom_init", "could not allocate blocks table" );
        free( result );
        fclose( f );
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    for(i = 0; i < dist_m; i++) {
        for(j = 0; j < dist_n; j++) {
            int rank;

            rc = fscanf( f, "%d", &rank );
            if ( rc < 1 ) {
                char message[300];
                snprintf(message, 300, "file '%s': could not read value at position (%d, %d)", dist_file, i, j );
                chameleon_error( "chameleon_getrankof_custom_init", message );
                free( result->blocks_dist );
                free( result );
                fclose( f );
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }

            if ( (rank < 0 ) || (rank >= np) )
            {
                char message[300];
                snprintf( message, 300, "file '%s': value %d at position (%d, %d) is invalid with %d processes",
                          dist_file, rank, i, j, np );
                chameleon_error( "chameleon_getrankof_custom_init", message );
                free( result->blocks_dist );
                free( result );
                fclose( f );
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }

            result->blocks_dist[j * dist_m + i] = rank;
        }
    }
    fclose(f);

    *custom_dist = result;
    return CHAMELEON_SUCCESS;
}

/**
 * @brief Destroys a custom distribution based on an external file.
 *
 * @param[int] dist
 *        The custom distribution to be destroyed
 *
 * @return CHAMELEON_SUCCESS on successful exit, CHAMELEON_ERR_UNALLOCATED otherwise.
 *
 */
int
chameleon_getrankof_custom_destroy( custom_dist_t **dist )
{
    if ((dist == NULL) || (*dist == NULL)) {
        chameleon_error("chameleon_getrankof_custom_destroy", "attempting to destroy a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    free((*dist)->blocks_dist);
    free(*dist);
    *dist = NULL;
    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to return MPI rank of block (m,n) in distribution
 * custom from dist file
 *
 * @param[in] desc matrix
 * @param[in] m row index of tile to consider
 * @param[in] n column index of tile to consider
 *
 * @return The rank of the tile at coordinate (m, n).
 */
int chameleon_getrankof_custom( const CHAM_desc_t *desc, int m, int n )
{
    custom_dist_t *dist = desc->get_rankof_init_arg;
    return dist->blocks_dist[(n % dist->dist_n) * dist->dist_m + (m % dist->dist_m)];
}

/**
 * @brief Return the address of the tile A( m, n ) in a tile storage.
 *
 * @WARNING The only mapping valid with this data storage is the 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The address of the tile A( m, n )
 *
 */
void* chameleon_getaddr_ccrb( const CHAM_desc_t *A, int m, int n )
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset  = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn ) );
    mm = mm / chameleon_desc_datadist_get_iparam(A, 0);
    nn = nn / chameleon_desc_datadist_get_iparam(A, 1);
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = (size_t)(A->bsiz) * (mm + (size_t)(A->llm1) * nn );
        else
            offset = A->A12 + ((size_t)(A->mb * (A->lln%A->nb)) * mm );
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((size_t)((A->llm%A->mb) * A->nb) * nn );
        else
            offset = A->A22;
    }

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 * @brief Return the address of the tile A( m, n ) in a column major storage.
 *
 * @WARNING The only mapping valid with this data storage is the 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The address of the tile A( m, n )
 *
 */
void *chameleon_getaddr_cm( const CHAM_desc_t *A, int m, int n )
{
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset  = chameleon_getaddr_cm_offset( A, m, n, A->llm );

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 * @brief Return the address of the tile A( m, m ) in a tile storage.
 *
 * @WARNING The only mapping valid with this data storage is the diagonal 2D block cyclic
 * storage as soon as multiple nodes are involved (see chameleon_getrankof_2d_diag()).
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row and column index of the tile.
 *
 * @param[in] n
 *        Unused.
 *
 * @return The address of the tile A( m, n )
 *
 */
void *chameleon_getaddr_diag( const CHAM_desc_t *A, int m, int n )
{
    (void)n;
    return chameleon_getaddr_ccrb( A, m, 0 );
}

/**
 * @brief Return the address of the tile A( m, m ) in a dynamic storage.
 *
 * Template that returns a NULL address for all tiles. This is used jointly with
 * the StarPU runtime system to allocate data on the fly and it can be used with
 * any data mapping.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile (unused).
 *
 * @param[in] m
 *        The row index of the tile (unused).
 *
 * @param[in] n
 *        The column index of the tile (unused).
 *
 * @return A null pointer
 *
 */
void *chameleon_getaddr_null( const CHAM_desc_t *A, int m, int n )
{
    (void)A; (void)m; (void)n;
    return NULL;
}

/**
 * @brief Return the leading dimension of the tile A( m, m ) stored in a tiled storage.
 *
 * This functions returns the leading dimension of the tile A( m, n ) and since
 * the tiled storage is compatible with any mapping strategy, this function is
 * also valid for any mapping.
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The leading dimension of the tile A( m, n ).
 *
 */
int chameleon_getblkldd_ccrb( const CHAM_desc_t *A, int m )
{
    int mm = m + A->i / A->mb;
    return ( ((mm+1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

/**
 * @brief Return the leading dimension of the tile A( m, m ) stored in a column major storage.
 *
 * This functions returns the leading dimension of the tile A( m, n ) gfor
 * column major storage. It can be only be used jointly with
 * chameleon_getaddr_cm() and chameleon_getrankof_2d().
 *
 * @param[in] A
 *        The matrix descriptor in which to find the tile.
 *
 * @param[in] m
 *        The row index of the tile.
 *
 * @param[in] n
 *        The column index of the tile.
 *
 * @return The leading dimension of the tile A( m, n ).
 *
 */
int chameleon_getblkldd_cm( const CHAM_desc_t *A, int m )
{
    (void)m;
    return A->llm;
}
