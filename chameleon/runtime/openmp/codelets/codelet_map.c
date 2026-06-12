/**
 *
 * @file openmp/codelet_map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map OpenMP codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2024-03-14
 *
 */
#include "chameleon_openmp.h"

typedef void (*omp_map_fct_t)( cham_uplo_t uplo, int m, int n, int ndata,
                               cham_map_data_t     *data,
                               cham_map_operator_t *op_fcts,
                               void                *op_args );

/*
 * @brief Open map function type with one tile
 */
static inline void omp_map_one_R( cham_uplo_t uplo, int m, int n, int ndata,
                                  cham_map_data_t     *data,
                                  cham_map_operator_t *op_fcts,
                                  void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );

    assert( ndata == 1 );
    assert( A->access == ChamR );

#pragma omp task depend( in:tileA[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 1, descA, tileA );
    }
}

static inline void omp_map_one_W( cham_uplo_t uplo, int m, int n, int ndata,
                                  cham_map_data_t     *data,
                                  cham_map_operator_t *op_fcts,
                                  void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );

    assert( ndata == 1 );
    assert( A->access == ChamW );

#pragma omp task depend( out:tileA[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 1, descA, tileA );
    }
}

static inline void omp_map_one_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                   cham_map_data_t     *data,
                                   cham_map_operator_t *op_fcts,
                                   void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );

    assert( ndata == 1 );
    assert( A->access == ChamRW );

#pragma omp task depend( in:tileA[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 1, descA, tileA );
    }
}

/*
 * @brief Open map function type with two tiles
 */
static inline void omp_map_two_R_R( cham_uplo_t uplo, int m, int n, int ndata,
                                    cham_map_data_t     *data,
                                    cham_map_operator_t *op_fcts,
                                    void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamR );
    assert( B->access == ChamR );

#pragma omp task depend( in:tileA[0] ) depend( in:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_R_W( cham_uplo_t uplo, int m, int n, int ndata,
                                    cham_map_data_t     *data,
                                    cham_map_operator_t *op_fcts,
                                    void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamR );
    assert( B->access == ChamW );

#pragma omp task depend( in:tileA[0] ) depend( out:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_R_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                     cham_map_data_t     *data,
                                     cham_map_operator_t *op_fcts,
                                     void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamR );
    assert( B->access == ChamRW );

#pragma omp task depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_W_R( cham_uplo_t uplo, int m, int n, int ndata,
                                    cham_map_data_t     *data,
                                    cham_map_operator_t *op_fcts,
                                    void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamW );
    assert( B->access == ChamR );

#pragma omp task depend( out:tileA[0] ) depend( in:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_W_W( cham_uplo_t uplo, int m, int n, int ndata,
                                    cham_map_data_t     *data,
                                    cham_map_operator_t *op_fcts,
                                    void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );

#pragma omp task depend( out:tileA[0] ) depend( out:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_W_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                     cham_map_data_t     *data,
                                     cham_map_operator_t *op_fcts,
                                     void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamW );
    assert( B->access == ChamRW );

#pragma omp task depend( out:tileA[0] ) depend( inout:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_RW_R( cham_uplo_t uplo, int m, int n, int ndata,
                                     cham_map_data_t     *data,
                                     cham_map_operator_t *op_fcts,
                                     void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamRW );
    assert( B->access == ChamR );

#pragma omp task depend( inout:tileA[0] ) depend( in:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_RW_W( cham_uplo_t uplo, int m, int n, int ndata,
                                     cham_map_data_t     *data,
                                     cham_map_operator_t *op_fcts,
                                     void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamRW );
    assert( B->access == ChamW );

#pragma omp task depend( inout:tileA[0] ) depend( out:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

static inline void omp_map_two_RW_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                      cham_map_data_t     *data,
                                      cham_map_operator_t *op_fcts,
                                      void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );

    assert( ndata == 2 );
    assert( A->access == ChamRW );
    assert( B->access == ChamRW );

#pragma omp task depend( inout:tileA[0] ) depend( inout:tileB[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 2, descA, tileA, descB, tileB );
    }
}

/*
 * @brief Open map function type with three tiles
 */
static inline void omp_map_three_R_R_R( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamR );
    assert( C->access == ChamR );

#pragma omp task depend( in:tileA[0] ) depend( in:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_R_W( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamW );

#pragma omp task depend( in:tileA[0] ) depend( in:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_R_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamR );
    assert( C->access == ChamRW );

#pragma omp task depend( in:tileA[0] ) depend( in:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_W_R( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( in:tileA[0] ) depend( out:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_W_W( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( in:tileA[0] ) depend( out:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_W_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamW );
    assert( C->access == ChamRW );

#pragma omp task depend( in:tileA[0] ) depend( out:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_RW_R( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamRW );
    assert( C->access == ChamR );

#pragma omp task depend( in:tileA[0] ) depend( inout:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_RW_W( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamRW );
    assert( C->access == ChamW );

#pragma omp task depend( in:tileA[0] ) depend( inout:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_R_RW_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamR );
    assert( B->access == ChamRW );
    assert( C->access == ChamRW );

#pragma omp task depend( in:tileA[0] ) depend( inout:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_R_R( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamR );
    assert( C->access == ChamR );

#pragma omp task depend( out:tileA[0] ) depend( in:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_R_W( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );

#pragma omp task depend( out:tileA[0] ) depend( in:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_R_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamR );
    assert( C->access == ChamRW );

#pragma omp task depend( out:tileA[0] ) depend( in:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_W_R( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( out:tileA[0] ) depend( out:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_W_W( cham_uplo_t uplo, int m, int n, int ndata,
                                        cham_map_data_t     *data,
                                        cham_map_operator_t *op_fcts,
                                        void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( out:tileA[0] ) depend( out:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_W_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamW );
    assert( C->access == ChamRW );

#pragma omp task depend( out:tileA[0] ) depend( out:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_RW_R( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamRW );
    assert( C->access == ChamR );

#pragma omp task depend( out:tileA[0] ) depend( inout:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_RW_W( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamRW );
    assert( C->access == ChamW );

#pragma omp task depend( out:tileA[0] ) depend( inout:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_W_RW_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamW );
    assert( B->access == ChamRW );
    assert( C->access == ChamRW );

#pragma omp task depend( out:tileA[0] ) depend( inout:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_R_R( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamR );
    assert( C->access == ChamR );

#pragma omp task depend( inout:tileA[0] ) depend( in:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_R_W( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamW );

#pragma omp task depend( inout:tileA[0] ) depend( in:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_R_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamR );
    assert( C->access == ChamRW );

#pragma omp task depend( inout:tileA[0] ) depend( in:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_W_R( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( inout:tileA[0] ) depend( out:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_W_W( cham_uplo_t uplo, int m, int n, int ndata,
                                         cham_map_data_t     *data,
                                         cham_map_operator_t *op_fcts,
                                         void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamW );
    assert( C->access == ChamR );

#pragma omp task depend( inout:tileA[0] ) depend( out:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_W_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamW );
    assert( C->access == ChamRW );

#pragma omp task depend( inout:tileA[0] ) depend( out:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_RW_R( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamRW );
    assert( C->access == ChamR );

#pragma omp task depend( inout:tileA[0] ) depend( inout:tileB[0] ) depend( in:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_RW_W( cham_uplo_t uplo, int m, int n, int ndata,
                                          cham_map_data_t     *data,
                                          cham_map_operator_t *op_fcts,
                                          void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamRW );
    assert( C->access == ChamW );

#pragma omp task depend( inout:tileA[0] ) depend( inout:tileB[0] ) depend( out:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static inline void omp_map_three_RW_RW_RW( cham_uplo_t uplo, int m, int n, int ndata,
                                           cham_map_data_t     *data,
                                           cham_map_operator_t *op_fcts,
                                           void                *op_args )
{
    cham_map_data_t   *A     = data;
    const CHAM_desc_t *descA = A->desc;
    CHAM_tile_t       *tileA = descA->get_blktile( descA, m, n );
    cham_map_data_t   *B     = data + 1;
    const CHAM_desc_t *descB = B->desc;
    CHAM_tile_t       *tileB = descB->get_blktile( descB, m, n );
    cham_map_data_t   *C     = data + 2;
    const CHAM_desc_t *descC = C->desc;
    CHAM_tile_t       *tileC = descC->get_blktile( descC, m, n );

    assert( ndata == 3 );
    assert( A->access == ChamRW );
    assert( B->access == ChamRW );
    assert( C->access == ChamRW );

#pragma omp task depend( inout:tileA[0] ) depend( inout:tileB[0] ) depend( inout:tileC[0] )
    {
        op_fcts->cpufunc( op_args, uplo, m, n, 3, descA, tileA, descB, tileB, descC, tileC );
    }
}

static omp_map_fct_t omp_map_one[3] = {
    omp_map_one_R, omp_map_one_W, omp_map_one_RW
};

static omp_map_fct_t omp_map_two[3][3] = {
    { omp_map_two_R_R,  omp_map_two_R_W,  omp_map_two_R_RW },
    { omp_map_two_W_R,  omp_map_two_W_W,  omp_map_two_W_RW },
    { omp_map_two_RW_R, omp_map_two_RW_W, omp_map_two_RW_RW }
};

static omp_map_fct_t omp_map_three[3][3][3] = {
    {   { omp_map_three_R_R_R,   omp_map_three_R_R_W,   omp_map_three_R_R_RW   },
        { omp_map_three_R_W_R,   omp_map_three_R_W_W,   omp_map_three_R_W_RW   },
        { omp_map_three_R_RW_R,  omp_map_three_R_RW_W,  omp_map_three_R_RW_RW  }  },
    {   { omp_map_three_W_R_R,   omp_map_three_W_R_W,   omp_map_three_W_R_RW   },
        { omp_map_three_W_W_R,   omp_map_three_W_W_W,   omp_map_three_W_W_RW   },
        { omp_map_three_W_RW_R,  omp_map_three_W_RW_W,  omp_map_three_W_RW_RW  }  },
    {   { omp_map_three_RW_R_R,  omp_map_three_RW_R_W,  omp_map_three_RW_R_RW  },
        { omp_map_three_RW_W_R,  omp_map_three_RW_W_W,  omp_map_three_RW_W_RW  },
        { omp_map_three_RW_RW_R, omp_map_three_RW_RW_W, omp_map_three_RW_RW_RW }  }
};


void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    switch( ndata ) {
    case 1:
        omp_map_one[ data[0].access - 1 ](
            uplo, m, n, ndata, data, op_fcts, op_args );
        break;

    case 2:
        omp_map_two[ data[0].access - 1 ][ data[1].access - 1 ](
            uplo, m, n, ndata, data, op_fcts, op_args );
        break;

    case 3:
        omp_map_three[ data[0].access - 1 ][ data[1].access - 1 ][ data[2].access - 1 ](
            uplo, m, n, ndata, data, op_fcts, op_args );
        break;
    }
}
