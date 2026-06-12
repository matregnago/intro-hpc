/**
 *
 * @file link_chameleon.c
 *
 * @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon link_chameleon example
 *
 * @version 1.4.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2025-12-19
 *
 */
#include <chameleon.h>
#include <stdlib.h>

/**
 * @brief Matrix multiplication example to check library linking with Chameleon
 */
int main() {
    /* number of threads */
    int NT = -1;
    /* number of gpus */
    int NG = 0;
    /* size of the matrices */
    int N = 16;
    /* size of the tiles */
    int NB = 2;
    /* seeds for A and B initialization */
    int seedA = 33;
    int seedB = 47;

    CHAMELEON_Init( NT, NG );
    CHAMELEON_Set( CHAMELEON_TILE_SIZE, NB );
    int NMPI = CHAMELEON_Comm_size();
    int GRID_P = NMPI;
    int GRID_Q = 1;

    CHAM_desc_t *descA = NULL, *descB = NULL, *descC = NULL;
    CHAMELEON_Desc_Create(&descA, NULL, ChamRealFloat,
        NB, NB, NB*NB, N, N, 0, 0, N, N, GRID_P, GRID_Q);
    CHAMELEON_Desc_Create(&descB, NULL, ChamRealFloat,
        NB, NB, NB*NB, N, N, 0, 0, N, N, GRID_P, GRID_Q);
    CHAMELEON_Desc_Create(&descC, NULL, ChamRealFloat,
        NB, NB, NB*NB, N, N, 0, 0, N, N, GRID_P, GRID_Q);

    CHAMELEON_splrnt_Tile( descA, 33 );
    CHAMELEON_splrnt_Tile( descB, 47 );
    CHAMELEON_sgemm_Tile(ChamNoTrans, ChamNoTrans, 1.0, descA, descB, 0.0, descC);

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    CHAMELEON_Finalize();

    return EXIT_SUCCESS;
}
