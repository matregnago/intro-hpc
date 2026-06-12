/**
 *
 * @file common.c
 *
 * @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 0.1.0
 * @author Mathieu Faverge
 * @date 2017-04-27
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "common.h"

#if defined(LIBHQR_HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(LIBHQR_HAVE_GETOPT_H) */

#ifndef max
#define max( _a, _b ) ((_a) < (_b) ? (_a) : (_b))
#endif

/**
 * @brief Print usage details of the testings
 */
static void
print_usage(void)
{
    fprintf(stderr,
            "Mandatory argument:\n"
            " -N                : Dimension (NT) of the matrices (required)\n"
            "Optional arguments:\n"
            " -P --grid-rows    : Rows (P) in the PxQ process grid (default: 1)\n"
            " -Q --grid-cols    : Cols (Q) in the PxQ process grid (default: 1)\n"
            " -c --cores        : Number of cores per node for automatically configured  trees (default: 1)\n"
            "\n"
            " -M                : Dimension (MT) of the matrices (default: NT)\n"
            "\n"
            " -a --qr_a         : Size of TS domain. (default: -1)\n"
            " -p --qr_p         : Size of the high level tree for distributed mode (default: -1)\n"
            " -d --domino       : Enable/Disable the domino between upper and lower trees. (default: -1)\n"
            " -r --tsrr         : Enable/Disable the round-robin on TS domain (under dev.). (default: Disabled)\n"
            " -l --treel        : Tree used for low level reduction inside nodes (default: -1).\n"
            " -L --treeh        : Tree used for high level reduction between nodes, only if qr_p > 1 (default: -1).\n"
            "                      (0: Flat, 1: Greedy, 2: Fibonacci, 3: Binary, 4: Replicated greedy)\n"
            "\n"
            " -x --check        : verify the results\n"
            " -X                  Check the results and print statistics for cdash\n"
            " -v --verbose      : extra verbose output\n"
            " -h --help         : this message\n"
            "\n"
            );
}

#define GETOPT_STRING "M:N:P:Q:c:a:p:drl:L:xXv::h"
#if defined(LIBHQR_HAVE_GETOPT_LONG)
static struct option long_options[] =
{
    /* Generic Options */
    {"M",           required_argument,  0, 'M'},
    {"N",           required_argument,  0, 'N'},
    {"grid-rows",   required_argument,  0, 'P'},
    {"grid-cols",   required_argument,  0, 'Q'},
    {"cores",       required_argument,  0, 'c'},

    /* HQR options */
    {"qr_a",        required_argument,  0, 'a'},
    {"qr_p",        required_argument,  0, 'p'},
    {"domino",      no_argument,        0, 'd'},
    {"tsrr",        no_argument,        0, 'r'},
    {"treel",       required_argument,  0, 'l'},
    {"treeh",       required_argument,  0, 'L'},

    /* Auxiliary options */
    {"check",       no_argument,        0, 'x'},
    {"X",           no_argument,        0, 'X'},
    {"verbose",     optional_argument,  0, 'v'},
    {"help",        no_argument,        0, 'h'},
    {0, 0, 0, 0}
};
#endif  /* defined(LIBHQR_HAVE_GETOPT_LONG) */

/**
 * @brief Set default values
 */
static void
iparam_default(int *iparam)
{
    /* Just in case someone forget to add the initialization :) */
    memset(iparam, 0, IPARAM_SIZEOF * sizeof(int));
    iparam[IPARAM_NNODES]       = 1;
    iparam[IPARAM_NCORES]       = 1;
    iparam[IPARAM_P]            = 1;
    iparam[IPARAM_Q]            = 1;
    iparam[IPARAM_MT]           = -'N';
    iparam[IPARAM_NT]           = 1;
    iparam[IPARAM_CHECK]        = 0;
    iparam[IPARAM_VERBOSE]      = 0;
    iparam[IPARAM_LOWLVL_TREE]  = -1;
    iparam[IPARAM_HIGHLVL_TREE] = -1;
    iparam[IPARAM_QR_TS_SZE]    = -1;
    iparam[IPARAM_QR_HLVL_SZE]  = -'P';
    iparam[IPARAM_QR_DOMINO]    = -1;
    iparam[IPARAM_QR_TSRR]      = 0;
}

/**
 * @brief Parse options
 */
void
parse_arguments(int *_argc, char ***_argv, int *iparam)
{
    int opt = 0;
    int c;
    int argc = *_argc;
    char **argv = *_argv;

    iparam_default( iparam );

    do {
#if defined(LIBHQR_HAVE_GETOPT_LONG)
        c = getopt_long_only(argc, argv, "",
                             long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(LIBHQR_HAVE_GETOPT_LONG) */

        switch(c)
        {
        case 'c':
            iparam[IPARAM_NCORES] = atoi(optarg);
            break;
        case 'P':
            iparam[IPARAM_P] = atoi(optarg);
            break;
        case 'Q':
            iparam[IPARAM_Q] = atoi(optarg);
            break;
        case 'M':
            iparam[IPARAM_MT] = atoi(optarg);
            break;
        case 'N':
            iparam[IPARAM_NT] = atoi(optarg);
            break;
        case 'x':
            iparam[IPARAM_CHECK] = 1;
            iparam[IPARAM_VERBOSE] = max(2, iparam[IPARAM_VERBOSE]);
            break;
        case 'X':
            iparam[IPARAM_CHECK] = 2;
            iparam[IPARAM_VERBOSE] = max(2, iparam[IPARAM_VERBOSE]);
            break;

            /* HQR parameters */
        case 'a':
            iparam[IPARAM_QR_TS_SZE] = atoi(optarg);
            break;
        case 'p':
            iparam[IPARAM_QR_HLVL_SZE] = atoi(optarg);
            break;
        case 'd':
            iparam[IPARAM_QR_DOMINO] = 1;
            break;
        case 'r':
            iparam[IPARAM_QR_TSRR] = 1;
            break;
        case 'l':
            iparam[IPARAM_LOWLVL_TREE]  = atoi(optarg);
            break;
        case 'L':
            iparam[IPARAM_HIGHLVL_TREE] = atoi(optarg);
            break;

        case 'v':
            if(optarg) {
                iparam[IPARAM_VERBOSE] = atoi(optarg);
            }
            else {
                iparam[IPARAM_VERBOSE] = 2;
            }
            break;

        case 'h':
            print_usage();
            exit(0);

        case '?': /* getopt_long already printed an error message. */
            exit(1);

        default:
            break; /* Assume anything else is parsec/mpi stuff */
        }
    } while(-1 != c);

    if(-'N' == iparam[IPARAM_MT]) {
        iparam[IPARAM_MT] = iparam[IPARAM_NT];
    }
    if(-'P' == iparam[IPARAM_QR_HLVL_SZE]) {
        iparam[IPARAM_QR_HLVL_SZE] = iparam[IPARAM_P];
    }
    if(-'Q' == iparam[IPARAM_QR_HLVL_SZE]) {
        iparam[IPARAM_QR_HLVL_SZE] = iparam[IPARAM_Q];
    }

    iparam[IPARAM_NNODES] = iparam[IPARAM_P] * iparam[IPARAM_Q];
}
