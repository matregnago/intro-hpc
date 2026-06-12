/**
 *
 * @file chameleon_tasks_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon complex mixed precision elementary tasks header
 *
 * @version 1.4.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions mixed zc -> ds
 *
 */
#ifndef _chameleon_tasks_zc_h_
#define _chameleon_tasks_zc_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );

#endif /* _chameleon_tasks_zc_h_ */
