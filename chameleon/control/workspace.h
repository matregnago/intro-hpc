/**
 *
 * @file workspace.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon workspace header
 *
 * @version 1.4.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_workspace_h_
#define _chameleon_workspace_h_

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
int chameleon_alloc_ibnb_tile(int M, int N, cham_tasktype_t func, int type, CHAM_desc_t **desc, int p, int q);
int chameleon_alloc_ipiv(int M, int N, cham_tasktype_t func, int type, CHAM_desc_t **desc, void **IPIV, int p, int q);

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_workspace_h_ */
