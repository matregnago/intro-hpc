/**
 *
 * @file primes.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon InPlaceTransformation prime numbers module header
 * @version 1.4.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2024-02-18
 *
 * This work is the implementation of an inplace transformation
 * based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 * and its fortran implementation.
 *
 */
#ifndef _primes_h_
#define _primes_h_

#define IMBALANCE_THRESHOLD 10
#define PWR_MAXSIZE   32
#define PRIME_MAXSIZE 10
#define SIZE_MG       1024
#define SIZE_LEADERS  1023

#ifndef min
#define min(a,b) ((a<b)?a:b)
#endif

#ifndef max
#define max(a,b) ((a>b)?a:b)
#endif

struct primedec
{
  int p;
  int e;
  int pe;
};

typedef struct primedec primedec_t;

int lcm(int a, int b);
int gcd(int a, int b);
int modpow(int x, int n, int m);
void factor(int n, primedec_t *pr, int *nf);

int     minloc(int n, int *T);
int64_t maxval(int n, int *T);
int64_t sum   (int n, int *T);

#endif /* _primes_h_ */
