/**
 *
 * @file starpu/codelet_gemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon gemm StarPU codelet
 *
 * @version 1.4.0
 * @author Mathieu Faverge
 * @date 2024-10-18
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelets.h"
#define PRECISION_z
#include "runtime_codelet_z.h"
#undef PRECISION_z
#define PRECISION_d
#include "runtime_codelet_d.h"
#undef PRECISION_d
#define PRECISION_c
#include "runtime_codelet_c.h"
#undef PRECISION_c
#define PRECISION_s
#include "runtime_codelet_s.h"
#undef PRECISION_s

void
INSERT_TASK_gemm( const RUNTIME_option_t *options,
                  cham_trans_t transA, cham_trans_t transB,
                  int m, int n, int k, int nb,
                  double alpha, const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *B, int Bm, int Bn,
                  double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    struct starpu_codelet *codelet = NULL;
    callback_fct_t         cl_gemm_callback = NULL;

    /* if ( alpha == 0. ) { */
    /*     INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb, */
    /*                          beta, C, Cm, Cn ); */
    /*     return; */
    /* } */

    void          *clargs = NULL;
    cham_access_t  cham_accessC;
    int            stpu_accessC;
    int            exec = 0;
    size_t         argssize = 0;
    const char    *cl_name = "Xgemm";
    CHAM_tile_t   *tileC;
    cham_flttype_t Cflttype;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Reduce the C access if needed */
    cham_accessC = ( beta == 0. ) ? ChamW    : ChamRW;
    stpu_accessC = ( beta == 0. ) ? STARPU_W : (STARPU_RW | ((beta == 1.) ? STARPU_COMMUTE : 0));

    tileC = C->get_blktile( C, Cm, Cn );
    Cflttype = tileC->flttype;

    switch( Cflttype ) {
#if defined(CHAMELEON_PREC_Z)
    case ChamComplexDouble:
        codelet  = &cl_zgemm;
        cl_gemm_callback = cl_zgemm_callback;
        if ( exec ) {
            struct cl_zgemm_args_s *cl_zargs;
            cl_zargs = malloc( sizeof( struct cl_zgemm_args_s ) );
            cl_zargs->transA = transA;
            cl_zargs->transB = transB;
            cl_zargs->m      = m;
            cl_zargs->n      = n;
            cl_zargs->k      = k;
            cl_zargs->alpha  = alpha;
            cl_zargs->beta   = beta;
            clargs = (void*)cl_zargs;
            argssize = sizeof( struct cl_zgemm_args_s );
        }
        break;
#endif
#if defined(CHAMELEON_PREC_C)
    case ChamComplexSingle:
        codelet  = &cl_cgemm;
        cl_gemm_callback = cl_cgemm_callback;
        if ( exec ) {
            struct cl_cgemm_args_s *cl_cargs;
            cl_cargs = malloc( sizeof( struct cl_cgemm_args_s ) );
            cl_cargs->transA = transA;
            cl_cargs->transB = transB;
            cl_cargs->m      = m;
            cl_cargs->n      = n;
            cl_cargs->k      = k;
            cl_cargs->alpha  = alpha;
            cl_cargs->beta   = beta;
            clargs = (void*)cl_cargs;
            argssize = sizeof( struct cl_cgemm_args_s );
        }
        break;
#endif
#if defined(CHAMELEON_PREC_D)
    case ChamRealDouble:
        codelet  = &cl_dgemm;
        cl_gemm_callback = cl_dgemm_callback;
        if ( exec ) {
            struct cl_dgemm_args_s *cl_dargs;
            cl_dargs = malloc( sizeof( struct cl_dgemm_args_s ) );
            cl_dargs->transA = transA;
            cl_dargs->transB = transB;
            cl_dargs->m      = m;
            cl_dargs->n      = n;
            cl_dargs->k      = k;
            cl_dargs->alpha  = alpha;
            cl_dargs->beta   = beta;
            clargs = (void*)cl_dargs;
            argssize = sizeof( struct cl_dgemm_args_s );
        }
        break;
#endif
#if defined(CHAMELEON_PREC_S)
    case ChamRealSingle:
        codelet  = &cl_sgemm;
        cl_gemm_callback = cl_sgemm_callback;
        if ( exec ) {
            struct cl_sgemm_args_s *cl_sargs;
            cl_sargs = malloc( sizeof( struct cl_sgemm_args_s ) );
            cl_sargs->transA = transA;
            cl_sargs->transB = transB;
            cl_sargs->m      = m;
            cl_sargs->n      = n;
            cl_sargs->k      = k;
            cl_sargs->alpha  = alpha;
            cl_sargs->beta   = beta;
            clargs = (void*)cl_sargs;
            argssize = sizeof( struct cl_sgemm_args_s );
        }
        break;
#endif
#if (defined(CHAMELEON_PREC_D) || defined(CHAMELEON_PREC_S)) && defined(CHAMELEON_USE_CUDA)
    case ChamRealHalf:
        codelet  = &cl_hgemm;
        cl_gemm_callback = cl_hgemm_callback;
        if ( exec ) {
            struct cl_hgemm_args_s *cl_hargs;
            cl_hargs = malloc( sizeof( struct cl_hgemm_args_s ) );
            cl_hargs->transA = transA;
            cl_hargs->transB = transB;
            cl_hargs->m      = m;
            cl_hargs->n      = n;
            cl_hargs->k      = k;
            cl_hargs->alpha  = alpha;
            cl_hargs->beta   = beta;
            clargs = (void*)cl_hargs;
            argssize = sizeof( struct cl_hgemm_args_s );
        }
        break;
#endif
    default:
        fprintf( stderr, "INSERT_TASK_gemm: Unknown datatype %d (Mixed=%3s, Type=%d, Size=%d\n",
                 Cflttype, cham_is_mixed(Cflttype) ? "Yes" : "No",
                 cham_get_ftype(Cflttype), cham_get_arith(Cflttype) );
        return;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 3,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        codelet,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, argssize,

        /* Task handles */
        STARPU_R,     RUNTIME_data_getaddr_withconversion( options, ChamR,        Cflttype, A, Am, An ),
        STARPU_R,     RUNTIME_data_getaddr_withconversion( options, ChamR,        Cflttype, B, Bm, Bn ),
        stpu_accessC, RUNTIME_data_getaddr_withconversion( options, cham_accessC, Cflttype, C, Cm, Cn ),

        /* Common task arguments */
        INSERT_TASK_COMMON_TASK_PARAMS( gemm ),
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
    return;
}
