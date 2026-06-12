/**
 *
 * @file gensvg.c
 *
 * File for algorithm of treewalking.
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                 Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#include "libhqr_internal.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH  50
#define HEIGHT 50
#define SIZE   100

/*
 * Global array for color
 */
char *colortree[] = {"red", "blue", "green", "orange", "cyan", "purple", "yellow" };
#define NBCOLORS (sizeof( colortree ) / sizeof( char* ))

/**
 * @brief Write the svg header
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 */
static void
drawsvg_header( FILE *file )
{
    int rc;

    rc = fprintf(file,
                 "<?xml version=\"1.0\" standalone=\"no\"?>\n"
                 "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
                 "<svg version=\"1.1\" \n xmlns=\"http://www.w3.org/2000/svg\" >\n");

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_header)\n");
    }
    return;
}

/**
 * @brief Write the box for the top node of a TS elimination
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] x
 *             The abscissa coordinate of the node
 * @param[in] y
 *             The ordinate coordinate of the node
 * @param[in] w
 *             The width of the node
 * @param[in] h
 *             The height of the node
 */
static void
drawsvg_top_TS( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf( file, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n",
                  x, y, w, h, colortree[k%NBCOLORS] );

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_top_TS)\n");
    }
    return;
}

/**
 * @brief Write the box for the bottom node of a TS elimination
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] x
 *             The abscissa coordinate of the node
 * @param[in] y
 *             The ordinate coordinate of the node
 * @param[in] w
 *             The width of the node
 * @param[in] h
 *             The height of the node
 */
static void
drawsvg_bot_TS( FILE *file, int k, int x, int y, int w, int h )
{
    int rc, x2, y2, w2, h2;

    rc = fprintf( file, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n",
                  x, y, w, h, colortree[k%NBCOLORS] );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TS)\n");
        return;
    }

    x2 = x + (w / 4);
    y2 = y + (h / 4);
    w2 = (w / 2);
    h2 = (h / 2);

    rc = fprintf( file, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill =\"white\"/> \n",
                  x2, y2, w2, h2 );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TS)\n");
    }
    return;
}

/**
 * @brief Write the circle for the top node of a TT elimination
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] x
 *             The abscissa coordinate of the node
 * @param[in] y
 *             The ordinate coordinate of the node
 * @param[in] w
 *             The width of the node
 * @param[in] h
 *             The height of the node
 */
static void
drawsvg_top_TT( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf( file, "<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" /> \n",
                  x + w / 2, y + h / 2, w / 2, colortree[k%NBCOLORS] );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_top_TT)\n");
    }
    return;
}

/**
 * @brief Write the circle for the bottom node of a TT elimination
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] x
 *             The abscissa coordinate of the node
 * @param[in] y
 *             The ordinate coordinate of the node
 * @param[in] w
 *             The width of the node
 * @param[in] h
 *             The height of the node
 */
static void
drawsvg_bot_TT( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf( file, "<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" /> \n",
                  x + w / 2, y + h / 2, w / 2, colortree[k%NBCOLORS] );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TT)\n");
        return;
    }

    rc = fprintf( file, "<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"white\" /> \n",
                  x + w / 2, y + h / 2, w / 4 );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TT)\n");
    }
    return;
}

/**
 * @brief Write an svg line from (x1,y1) to (x2,y2)
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] x1
 *             The abscissa coordinate of the origin node
 * @param[in] y1
 *             The ordinate coordinate of the origin node
 * @param[in] x2
 *             The abscissa coordinate of the destination node
 * @param[in] y2
 *             The ordinate coordinate of the destination node
 */
static void
drawsvg_line( FILE *file, int k, int x1, int y1, int x2, int y2 )
{
    int rc;
    rc = fprintf( file, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"fill:none;stroke:%s;stroke-width:2px;\"/>\n",
                  x1, y1, x2, y2, colortree[k%NBCOLORS] );

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_line)\n");
    }
    return;
}

/**
 * @brief Write the svg footer
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 */
static void
drawsvg_footer( FILE *file )
{
    int rc;
    rc = fprintf(file, "</svg>\n");
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_footer)\n");
    }
    return;
}

/**
 * @brief Draw the line corresponding to one TT/TS kernel
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] p
 *             The index of the pivot
 * @param[in] m
 *             The index of the tile being killed
 * @param[in] beg_p
 *             The last modification date of p
 * @param[in] beg_m
 *             The last modification date of m
 * @param[in] end
 *             The computed end of the operation
 */
static void
drawsvg_lines_rowm( FILE *file, int k,
                    int p, int m, int beg_p, int beg_m, int end )
{
    int yp, ym;
    int x, xp, xm;

    /* Row of the tiles */
    ym = SIZE + SIZE * m;
    yp = SIZE + SIZE * p;

    /* Starting position of the tiles */
    xm = (SIZE + (SIZE / 4)) + SIZE * beg_m;
    xp = (SIZE + (SIZE / 4)) + SIZE * beg_p;

    /* Final position of the tiles */
    x = SIZE + SIZE * end;

    /* Horizontal lines */
    drawsvg_line( file, k, xm, ym, x + (SIZE / 4), ym );
    drawsvg_line( file, k, xp, yp, x + (SIZE / 4), yp );

    /* Vertical line */
    drawsvg_line( file, k, x, ym, x, yp );
}

/**
 * @brief Draw all the lines related to a single step of the factorization
 * @param[in] qrtree
 *             The reduction tree structure
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in,out] tiles
 *             A temporary buffer of size qrtree->mt to store the list of tiles
 *             sorted by order of reduction
 * @param[in,out] step
 *             A buffer of size qrtree->mt with the date of last modification
 *             for each tile. On exit, the dates are updated with the performed
 *             updates.
 */
static void
drawsvg_lines_stepk( const libhqr_tree_t *qrtree, FILE *file,
                     int k, int *tiles, int *step )
{
    int i;
    int imax;

    /* Get order for step k */
    imax = libhqr_walk_stepk( qrtree, k, tiles );

    for(i = 0; i < imax; i++){
        int m = tiles[i];
        int p = qrtree->currpiv(qrtree, k, m);
        int end = libhqr_imax( step[p], step[m] ) + 1;

        /* Draw horizontal lines for rows p and m */
        drawsvg_lines_rowm( file, k, p, m, step[p], step[m], end );

        /* Update last time the rows p and m have been modified for future lines */
        step[m] = end;
        step[p] = end;
    }
}

/**
 * @brief Draw the two nodes related to one single reduction
 * @param[in,out] file
 *             The opened file decriptor on which the data is written
 * @param[in] k
 *             The factorization step that defines the color.
 * @param[in] type
 *             The type of reduction performed in the libhqr_type_e
 * @param[in] p
 *             The index of the pivot
 * @param[in] m
 *             The index of the tile being killed
 * @param[in] step_m
 *             The date at which the reduction occurs
 */
static void
drawsvg_nodes_rowm( FILE *file, int k,
                    libhqr_type_e type, int p, int m, int step_m )
{
    int x, yp, ym;
    x  = ((SIZE * 3) / 4) + SIZE * step_m;
    ym = ((SIZE * 3) / 4) + SIZE * m;
    yp = ((SIZE * 3) / 4) + SIZE * p;

    if ( type == LIBHQR_KILLED_BY_TS ) {
        drawsvg_top_TS(file, k, x, yp, WIDTH, HEIGHT );
        drawsvg_bot_TS(file, k, x, ym, WIDTH, HEIGHT );
    }
    else {
        drawsvg_top_TT(file, k, x, yp, WIDTH, HEIGHT );
        drawsvg_bot_TT(file, k, x, ym, WIDTH, HEIGHT );
    }
}

/**
 * @brief Draw a given reduction tree in an output svg file
 * @param[in] qrtree
 *             The reduction tree to draw
 * @param[in] filename
 *             The output filename of the svg
 */
void
libhqr_print_svg( const libhqr_tree_t *qrtree,
                  const char          *filename )
{
    FILE *file;
    int  *tiles, *steps;
    int   k, i;
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );
    steps = (int*)calloc( qrtree->mt, sizeof(int) );

    file = fopen(filename,"w+");

    drawsvg_header(file);
    for (k = 0; k < qrtree->nt; k++) {
        /* Drawing the lines */
        drawsvg_lines_stepk( qrtree, file, k, tiles, steps );

        /* Drawing the rectangles */
        for (i = 0; i < qrtree->mt; i++) {
            if (qrtree->currpiv( qrtree, k, i) < 0) {
                continue;
            }
            drawsvg_nodes_rowm( file, k,
                                qrtree->gettype(qrtree, k, i),
                                qrtree->currpiv(qrtree, k, i), i, steps[i] );
        }
    }
    drawsvg_footer(file);
    fclose(file);

    free(tiles);
    free(steps);
}
