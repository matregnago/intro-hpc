/**
 *
 * @file queue.c
 *
 * Queue module for the treewalk algorithm.
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                 Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-04-05
 *
 */
#include "libhqr_queue.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

libhqr_queue_tile_t *
libhqr_queue_tile_new(void)
{
    return NULL;
}

void
libhqr_queue_tile_push( libhqr_queue_tile_t **queue_tile, int numero )
{
    if (queue_tile != NULL)
    {
        libhqr_queue_tile_t *p_l;
        libhqr_queue_tile_t *p_p;

        assert( (*queue_tile) == NULL || (*queue_tile)->prev == NULL );
        p_l = *queue_tile;
        p_p = malloc (sizeof (*p_p));
        if (p_p != NULL)
        {
            p_p->numero = numero;
            p_p->next = p_l;
            p_p->prev = NULL;
            if (p_l != NULL) {
                p_l->prev = p_p;
            }
            *queue_tile = p_p;
        }
        else
        {
            fprintf (stderr, "Memoire insuffisante\n");
            exit (EXIT_FAILURE);
        }
    }
    return;
}

int
libhqr_queue_tile_head( libhqr_queue_tile_t **queue_tile )
{
    if ( ( queue_tile != NULL) &&
         (*queue_tile != NULL) )
    {
        return (*queue_tile)->numero;
    }
    return -1;
}

int
libhqr_queue_tile_pop ( libhqr_queue_tile_t **queue_tile )
{
    int ret = -1;
    if (queue_tile != NULL)
    {
        libhqr_queue_tile_t *p_l;
        libhqr_queue_tile_t *p_p;

        p_l = *queue_tile;
        p_p = p_l->next;
        if (p_p != NULL) {
            p_p->prev = NULL;
        }
        ret = p_l->numero;
        free (p_l);
        *queue_tile = p_p;
    }
    return ret;
}
