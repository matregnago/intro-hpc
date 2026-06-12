/**
 *
 * @file starpu/runtime_zprofiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU CHAMELEON_Complex64_t kernel progiling
 *
 * @version 1.4.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-09-17
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

void RUNTIME_zdisplay_allprofile()
{

    /* BLAS 3 */
    profiling_display_zgemm_info();
#if defined(PRECISION_z) || defined(PRECISION_c)
    profiling_display_zhemm_info();
    profiling_display_zher2k_info();
    profiling_display_zherk_info();
    profiling_display_zsytrf_nopiv_info();
#endif
    profiling_display_zsymm_info();
    profiling_display_zsyr2k_info();
    profiling_display_zsyrk_info();
    profiling_display_ztrmm_info();
    profiling_display_ztrsm_info();

    /* Lapack */
    profiling_display_zgelqt_info();
    profiling_display_zgeqrt_info();
    profiling_display_zgessm_info();
    profiling_display_zgetrf_incpiv_info();
    profiling_display_zgetrf_info();
    profiling_display_zgetrf_nopiv_info();
    profiling_display_zlauum_info();
    profiling_display_zpotrf_info();
    profiling_display_zssssm_info();
    profiling_display_ztplqt_info();
    profiling_display_ztpmlqt_info();
    profiling_display_ztpmqrt_info();
    profiling_display_ztpqrt_info();
    profiling_display_ztrtri_info();
    profiling_display_ztstrf_info();
    profiling_display_zunmlq_info();
    profiling_display_zunmqr_info();

    profiling_display_zlange_info();
    profiling_display_zlatro_info();
}

void RUNTIME_zdisplay_oneprofile( cham_tasktype_t kernel )
{
    switch( kernel ) {
        /* Blas 3 */
    case TASK_GEMM:         profiling_display_zgemm_info();         break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case TASK_HEMM:         profiling_display_zhemm_info();         break;
    case TASK_HER2K:        profiling_display_zher2k_info();        break;
    case TASK_HERK:         profiling_display_zherk_info();         break;
    case TASK_SYTRF_NOPIV:  profiling_display_zsytrf_nopiv_info();  break;
#endif
    case TASK_SYMM:         profiling_display_zsymm_info();         break;
    case TASK_SYR2K:        profiling_display_zsyr2k_info();        break;
    case TASK_SYRK:         profiling_display_zsyrk_info();         break;
    case TASK_TRMM:         profiling_display_ztrmm_info();         break;
    case TASK_TRSM:         profiling_display_ztrsm_info();         break;

        /* Lapack */
    case TASK_GELQT:        profiling_display_zgelqt_info();        break;
    case TASK_GEQRT:        profiling_display_zgeqrt_info();        break;
    case TASK_GESSM:        profiling_display_zgessm_info();        break;
    case TASK_GETRF:        profiling_display_zgetrf_info();        break;
    case TASK_GETRF_INCPIV: profiling_display_zgetrf_incpiv_info(); break;
    case TASK_GETRF_NOPIV:  profiling_display_zgetrf_nopiv_info();  break;
    case TASK_LAUUM:        profiling_display_zlauum_info();        break;
    case TASK_POTRF:        profiling_display_zpotrf_info();        break;
    case TASK_SSSSM:        profiling_display_zssssm_info();        break;
    case TASK_TRTRI:        profiling_display_ztrtri_info();        break;
    case TASK_TSTRF:        profiling_display_ztstrf_info();        break;

    case TASK_TPLQT:        profiling_display_ztplqt_info();        break;
    case TASK_TPMLQT:       profiling_display_ztpmlqt_info();        break;
    case TASK_TPMQRT:       profiling_display_ztpmqrt_info();        break;
    case TASK_TPQRT:        profiling_display_ztpqrt_info();        break;

    case TASK_UNMLQ:        profiling_display_zunmlq_info();        break;
    case TASK_UNMQR:        profiling_display_zunmqr_info();        break;

    case TASK_LANGE:        profiling_display_zlange_info();        break;
    case TASK_LATRO:        profiling_display_zlatro_info();        break;

    default:
        return;
    }
}

