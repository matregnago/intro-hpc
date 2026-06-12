/**
 *
 * @file libhqr_queue.h
 *
 * Queue module for the treewalk algorithm.
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#ifndef _libhqr_queue_h_
#define _libhqr_queue_h_

typedef struct libhqr_queue_tile_s {
    struct libhqr_queue_tile_s *prev;
    struct libhqr_queue_tile_s *next;
    int numero;
} libhqr_queue_tile_t;

libhqr_queue_tile_t *libhqr_queue_tile_new (void);
void libhqr_queue_tile_push  (libhqr_queue_tile_t **queue_tile, int numero);
int  libhqr_queue_tile_head  (libhqr_queue_tile_t **queue_tile);
int  libhqr_queue_tile_pop   (libhqr_queue_tile_t **queue_tile);
void libhqr_queue_tile_delete(libhqr_queue_tile_t **queue_tile);

#endif /* _libhqr_queue_h_ */
