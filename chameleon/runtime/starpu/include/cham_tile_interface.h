/**
 *
 * @file starpu/cham_tile_interface.h
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Header to describe the Chameleon tile interface in StarPU
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Gwenole Lucas
 * @author Ana Hourcau
 * @date 2024-07-17
 *
 */
#ifndef _cham_tile_interface_h_
#define _cham_tile_interface_h_

extern struct starpu_data_interface_ops starpu_interface_cham_tile_ops;
#define STARPU_CHAM_TILE_INTERFACE_ID starpu_interface_cham_tile_ops.interfaceid

struct starpu_cham_tile_interface_s;
typedef struct starpu_cham_tile_interface_s starpu_cham_tile_interface_t;

/**
 * Chameleon tile interface
 */
struct starpu_cham_tile_interface_s
{
    enum starpu_data_interface_id id; /**< Identifier of the interface           */
    uintptr_t      dev_handle;        /**< device handle of the matrix           */
    cham_flttype_t flttype;           /**< Type of the elements of the matrix    */
    size_t         allocsize;         /**< size actually currently allocated     */
    size_t         tilesize;          /**< size of the elements of the matrix    */
    CHAM_tile_t    tile;              /**< Internal tile structure used to store
                                           information on non memory home_node   */
};

void starpu_cham_tile_register( starpu_data_handle_t *handleptr,
                                int                   home_node,
                                CHAM_tile_t          *tile,
                                cham_flttype_t        flttype );

int    cti_handle_get_m        ( starpu_data_handle_t handle );
int    cti_handle_get_n        ( starpu_data_handle_t handle );
size_t cti_handle_get_allocsize( starpu_data_handle_t handle );

static inline CHAM_tile_t *
cti_interface_get( starpu_cham_tile_interface_t *interface )
{
    return &(interface->tile);
}

static inline CHAM_tile_t *
cti_handle_get( starpu_data_handle_t handle )
{
    starpu_cham_tile_interface_t *cham_tile_interface = (starpu_cham_tile_interface_t *)
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( cham_tile_interface->id == STARPU_CHAM_TILE_INTERFACE_ID,
                       "Error. The given data is not a cham_tile." );
#endif

    return &(cham_tile_interface->tile);
}

void starpu_cham_tile_interface_init();
void starpu_cham_tile_interface_fini();

#endif /* _cham_tile_interface_h_ */
