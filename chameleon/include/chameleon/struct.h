/**
 *
 * @file chameleon_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon structures
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @author Brieuc Nicolas
 * @date 2025-12-19
 *
 */
#ifndef _chameleon_struct_h_
#define _chameleon_struct_h_

#include "chameleon/config.h"
#include "chameleon/types.h"
#include "chameleon/constants.h"

BEGIN_C_DECLS

#define CHAMELEON_TILE_FULLRANK (1 << 0)
#define CHAMELEON_TILE_DESC     (1 << 1)
#define CHAMELEON_TILE_HMAT     (1 << 2)

#define CHAMELEON_MAX_DIMENSION 10

/**
 * @brief CHAMELEON structure to hold pivot informations for the LU factorization with partial pivoting
 */
typedef struct chameleon_pivot_s {
    int   blkm0;   /**> The row index of the first row in the tile where the pivot has been selected */
    int   blkidx;  /**> The relative row index in the tile where the pivot has been selected         */
    void *pivrow;  /**> The copy of the row with the selected pivot                                  */
    void *diagrow; /**> The copy of the diagonal row to permute                                      */
} CHAM_pivot_t;

typedef struct chameleon_tile_s {
#if defined(CHAMELEON_KERNELS_TRACE)
    char          *name;
#endif
    void          *mat;
    int            rank, m, n, ld;
    int8_t         format;
    cham_flttype_t flttype;
} CHAM_tile_t;

/**
 * @brief CHAMELEON interface representing the cpui workspace used for the permutation of LASWP and GETRF
 */
typedef struct chameleon_laswp_ws_s {
    int  *index;  /**< The array of indexes where the corresponding rows should be copied */
    void *rows;   /**< The array of rows to be moved depending of index                   */
    int   nindex; /**< The number of rows stored in the workspace                         */
    int   offset; /**< Offset to align the data after the indices array                   */
} CHAM_laswpws_t;

/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
struct chameleon_desc_s;
typedef struct chameleon_desc_s CHAM_desc_t;

typedef void*        (*blkaddr_fct_t)        ( const CHAM_desc_t*, int, int );
typedef int          (*blkldd_fct_t)         ( const CHAM_desc_t*, int );
typedef int          (*blkdim_fct_t)         ( const CHAM_desc_t*, int, cham_dim_t, int );
typedef int          (*blkrankof_fct_t)      ( const CHAM_desc_t*, int, int );
typedef int          (*datadist_access_fct_t)( const CHAM_desc_t*, int, ... );
typedef CHAM_tile_t* (*blktile_fct_t)        ( const CHAM_desc_t*, int, int );

/**
 * Data distribution type and acces functions
 */
/**
 * @brief Function discribing the indexed access to a 2D block cyclic data
 *        distribution holding the the grid size in the data distribution
 *        array [p,q] :
 *        0 returns p, 1 returns q.
 */
int chameleon_get_2d_block_cyclic(const CHAM_desc_t *desc, int i );

/**
 * @brief Function call which forward the index i to the access function
 */
int chameleon_desc_datadist_get_iparam( const CHAM_desc_t *desc, int i );

/**
 * @brief Data distribution type
 */
typedef struct cham_data_dist_s {
    datadist_access_fct_t get_distrib;    /**> function describing how to index the distribution array */
    int distrib_array_size;               /**> number of parameters stored in the distribution array   */
    int distrib[CHAMELEON_MAX_DIMENSION]; /**> array holding the parameters                            */
} cham_data_dist_t;

/**
 * @brief Function for initialising the data distribution
 */
void chameleon_desc_set_datadist( CHAM_desc_t *to, cham_data_dist_t *from );

struct chameleon_desc_s {
    char           *name; /**> Name of the descriptor for debug purpose */

    blktile_fct_t   get_blktile;     /**> function to get chameleon tiles address                     */
    blkaddr_fct_t   get_blkaddr;     /**> function to get chameleon tiles address                     */
    blkldd_fct_t    get_blkldd;      /**> function to get chameleon tiles leading dimension           */
    blkdim_fct_t    get_blkdim;      /**> function to get chameleon tiles dimension within algorithms */
    blkrankof_fct_t get_rankof;      /**> function to get chameleon tiles MPI rank                    */

    blkrankof_fct_t get_rankof_init;     /**> Pointer to the function used to set the data distribution
                                              (used only for checking purpose)                           */
    void           *get_rankof_init_arg; /**> Pointer to the data used by the data distribution function
                                              (used only to give the same prototype as get_blktile
                                              (should be modified) */

    CHAM_tile_t *tiles;  /**> pointer to the array of tiles descriptors  */
    void *mat;           /**> pointer to the beginning of the matrix     */
    size_t A21;          /**> pointer to the beginning of the matrix A21 */
    size_t A12;          /**> pointer to the beginning of the matrix A12 */
    size_t A22;          /**> pointer to the beginning of the matrix A22 */
    cham_storage_t styp; /**> storage layout of the matrix               */
    cham_flttype_t dtyp; /**> precision of the matrix                    */
    int mb;              /**> number of rows in a tile                   */
    int nb;              /**> number of columns in a tile                */
    int bsiz;            /**> size in elements including padding         */

    /* Matrix sizes in single rows/columns for the full problem */
    int i;            /**> row index to the beginning of the submatrix    */
    int j;            /**> column index to the beginning of the submatrix */
    int m;            /**> number of rows of the submatrix                */
    int n;            /**> number of columns of the submatrix             */
    int lm;           /**> number of rows of the entire matrix            */
    int ln;           /**> number of columns of the entire matrix         */

    /* Number of rows/columns tiles for the full problem */
    int mt;           /**> number of tile rows    of the submatrix - derived parameter     */
    int nt;           /**> number of tile columns of the submatrix - derived parameter     */
    int lmt;          /**> number of tile rows    of the entire matrix - derived parameter */
    int lnt;          /**> number of tile columns of the entire matrix - derived parameter */

    /* Distributed case */
    cham_data_dist_t *data_dist; /**> data distribution type used to retrieve the distributed layout      */
    int llm;                     /**> local number of rows         of the full matrix - derived parameter */
    int lln;                     /**> local number of columns      of the full matrix - derived parameter */
    int llm1;                    /**> local number of tile rows    of the A11  matrix - derived parameter */
    int lln1;                    /**> local number of tile columns of the A11  matrix - derived parameter */
    int llmt;                    /**> local number of tile rows    of the full matrix - derived parameter */
    int llnt;                    /**> local number of tile columns of the full matrix - derived parameter */

    int id;           /**> identification number of the descriptor                            */
    int occurences;   /**> identify main matrix desc (occurences=1) or                        */
                      /**> submatrix desc (occurences>1) to avoid unregistering               */
                      /**> GPU data twice                                                     */
    int sync;         /**> 1 if we the descriptor needs to be synced during flush - else 0    */
    int use_mat;      /**> 1 if we have a pointer to the overall data mat - else 0            */
    int alloc_mat;    /**> 1 if we handle the allocation of mat - else 0                      */
    int register_mat; /**> 1 if we have to register mat - else 0 (handled by the application) */
    int myrank;       /**> MPI rank of the descriptor                                         */
    int ooc;          /**> 1 if the matrix is not to fit in memory                            */
    int64_t mpitag;   /**> First MPI tag used by the descriptor                               */
    void *schedopt;   /**> scheduler (QUARK|StarPU) specific structure                        */
};

typedef struct chameleon_ipiv_s CHAM_ipiv_t;

typedef int (*blkdim_ipiv_fct_t)    ( const CHAM_ipiv_t*, int );
typedef int (*blkrankof_ipiv_fct_t) ( const CHAM_ipiv_t*, int, int );

/**
 * Bitwise enum to manage the three vectors stored in the CHAM_ipiv_t structure
 */
typedef enum cham_ipiv_ {
    CHAMIPIV_IPIV = (1 << 0),
    CHAMIPIV_PERM = (1 << 1),
    CHAMIPIV_INVP = (1 << 2),
} CHAM_ipiv_e;

/**
 *  CHAMELEON structure to hold ipiv/permutation informations for operation such as lapmt or laswp
 */
struct chameleon_ipiv_s {
    blkdim_ipiv_fct_t    get_blkdim;  /**> function to get chameleon tiles dimension within algorithms                                 */
    blkrankof_ipiv_fct_t get_rankof;  /**> function to get chameleon tiles MPI rank                                                    */
    int                 *data;        /**> Pointer to the data                                                                         */
    void                *ipiv;        /**> Opaque array of pointers for the runtimes to handle the ipiv array                          */
    void                *perm;        /**> Opaque array of pointers for the runtimes to handle the temporary permutation array         */
    void                *invp;        /**> Opaque array of pointers for the runtimes to handle the temporary inverse permutation array */
    int64_t              mpitag_ipiv; /**> Initial mpi tag values for the ipiv handles                                                 */
    int64_t              mpitag_perm; /**> Initial mpi tag values for the nextpiv handles                                              */
    int64_t              mpitag_invp; /**> Initial mpi tag values for the prevpiv handles                                              */
    int                  myrank;      /**> MPI rank of the descriptor                                                                  */
    int                  i;           /**> row index to the beginning of the submatrix                                                 */
    int                  m;           /**> The number of row in the vector ipiv                                                        */
    int                  mb;          /**> The number of row per block                                                                 */
    int                  mt;          /**> The number of tiles                                                                         */
    int                  withidx;     /**> Boolean to enable/disable the use of the {perm/invp}_idx instead of the classic ones        */
    int                  max_m;       /**> The number of rows/columns of the matrix to permute                                         */
    int                  max_mt;      /**> The number of row/clumn tiles of the matrix to permute                                      */
    int                  P;           /**> The number of processes per column on a tiled matrix                                        */
    int                  NP;          /**> The total number of processes                                                               */
};

/**
 * @brief CHAMELEON workspace structure used to apply permutations in kernels such as LASWP or GETRF
 */
typedef struct chameleon_perm_s {
    void           *ws;        /**< Workspace handle used for the row/column permutation */
    cham_side_t     side;      /**< Specifies the side of the permutation                */
    int64_t         mpitag_ws; /**< Initial MPI tag for the workspace handle             */
    int             NP;        /**< Total number of processes                            */
    int             m;         /**< Total number of rows of the matrix                   */
    int             n;         /**< Total number of columns of the matrix                */
    int             mb;        /**< Number of rows per tile                              */
    int             nb;        /**< Number of columns per tile                           */
    int             mt;        /**< Number of rows tile                                  */
    int             nt;        /**< Number of columns tile                               */
    cham_flttype_t  dtyp;      /**> Arithmetic used to store the rows/columns to swap    */
} CHAM_perm_t;

/**
 * Chameleon structure to hold pivot informations for the LU factorization with partial pivoting within a full panel
 * It includes two CHAM_pivot_t structures to manage the data dependencies between the columns.
 */
typedef struct chameleon_desc_pivot_s {
    void          *nextpiv;        /**> Opaque array of pointers for the runtimes to handle the pivot computation structure */
    void          *prevpiv;        /**> Opaque array of pointers for the runtimes to handle the pivot computation structure */
    int64_t        mpitag_nextpiv; /**> Initial mpi tag values for the nextpiv handles                                      */
    int64_t        mpitag_prevpiv; /**> Initial mpi tag values for the prevpiv handles                                      */
    int            P;              /**> The number of processes per column of the tiled matrix                              */
    int            Q;              /**> The number of processes per line of the tiled matrix                                */
    int            nb;             /**> The number of row per block                                                         */
    int            n;              /**> The number of column considered (must be updated for each panel) */
    cham_flttype_t dtyp;           /**> Arithmetic used to store the rows/columns to swap                                   */
} CHAM_desc_pivot_t;

static inline void *
CHAM_tile_get_ptr( const CHAM_tile_t *tile )
{
    if ( tile->format & CHAMELEON_TILE_DESC ) {
        return ((CHAM_desc_t*)(tile->mat))->mat;
    }
    return tile->mat;
}

static inline const char *
CHAM_tile_get_typestr( const CHAM_tile_t *tile )
{
    if ( tile->format & CHAMELEON_TILE_DESC ) {
        return "Desc";
    }

    if ( tile->format & CHAMELEON_TILE_HMAT ) {
        return "HMat";
    }

    return "Full";
}

END_C_DECLS

#endif /* _chameleon_struct_h_ */
