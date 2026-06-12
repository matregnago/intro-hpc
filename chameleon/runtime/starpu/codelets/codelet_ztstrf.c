/**
 *
 * @file starpu/codelet_ztstrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf StarPU codelet
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztstrf_cpu_func(void *descr[], void *cl_arg)
{
    CHAMELEON_starpu_ws_t *d_work;
    int m;
    int n;
    int ib;
    int nb;
    CHAM_tile_t *tileU;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileL;
    int *IPIV;
    CHAM_tile_t *tileW;
    int ldW;
    cham_bool_t check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    tileU = cti_interface_get(descr[0]);
    tileA = cti_interface_get(descr[1]);
    tileL = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]);

    starpu_codelet_unpack_args( cl_arg, &m, &n, &ib, &nb,
                                &IPIV, &d_work, &ldW, &check_info, &iinfo,
                                &sequence, &request );

    TCORE_ztstrf(m, n, ib, nb, tileU, tileA, tileL, IPIV, tileW->mat, ldW, &info);

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztstrf, cl_ztstrf_cpu_func)

void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo )
{
    (void)nb;
    CHAMELEON_starpu_ws_t *d_work = (CHAMELEON_starpu_ws_t*)(options->ws_host);

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(U, Um, Un);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_W(L, Lm, Ln);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        &cl_ztstrf,
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_VALUE,    &nb,                        sizeof(int),
        STARPU_RW,        RTBLKADDR(U, ChamComplexDouble, Um, Un),
        STARPU_RW,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,         RTBLKADDR(L, ChamComplexDouble, Lm, Ln),
        STARPU_VALUE,    &IPIV,                      sizeof(int*),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &d_work,                    sizeof(CHAMELEON_starpu_ws_t *),
        STARPU_VALUE,    &nb,                        sizeof(int),
        STARPU_VALUE,    &check_info,                sizeof(cham_bool_t),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( ztstrf ),
        0 );
}
