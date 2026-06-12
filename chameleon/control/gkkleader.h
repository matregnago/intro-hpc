/**
 *
 * @file gkkleader.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon InPlaceTransformation main module header
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 *
 * This work is the implementation of an inplace transformation
 * based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 * and its fortran implementation.
 *
 */
#ifndef _gkkleader_h_
#define _gkkleader_h_

int  GKK_doublingtable(int x, int m, int emax, int *dt);
int  GKK_modpow(int *dt, int e, int m);
int  GKK_primroot(int p, int e, primedec_t *pr_p1, int t_p1);
int  GKK_multorder(int n, int p, int e, int pe, primedec_t *pr_p1, int t_p1);
void GKK_prepare(int q, int n, primedec_t *pr_q, int *t,
                 primedec_t **pr_pi1, int *t_pi1, int *gi,
                 int *Nif, int *Kif, int *dif);
void GKK_L(int t, primedec_t *pr_q, int *fi, int *Kif, int *dif,
           int *Li, int *diLi, int *cl, int *nl);
void GKK_precompute_terms(int q, primedec_t *pr_q, int t, int *gi,
                          int *diLi, int *rp, int *Mg, int nMg);
void GKK_BalanceLoad(int thrdnbr, int *Tp, int *leaders, int nleaders, int L);
void GKK_output_tables(int m, int n, int q, primedec_t *pr_q, int t,
                       int *gi, int *Nif, int *Kif, int *dif);

int  GKK_getLeaderNbr(int me, int ne, int *nleaders, int **leaders);

#endif /* _gkkleader_h_ */
